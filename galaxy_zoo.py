import warnings
import numpy as np
import skimage as sk
import petrofit as pf
from pathlib import Path
from scipy import ndimage as ndi
from matplotlib import pyplot as plt


# Normalizes an image to avoid clippinggalaxy_zoo.py
def renorm(image):
    return np.clip(image, 0.0, 1.0)


# Compute the log scale of frequency after minimizing peaks for visualization
def log_scale(frequency_array):
    F = np.log(1 + np.abs(frequency_array))
    return renorm(F / F.max())


# Helper function that returns the FFT, including mag and angle
def fft_pair(image):
    F = np.fft.fft2(image)
    return F, np.abs(F), np.angle(F)


# takes 2d array of image, returns after adjusting gamma and cleaning sky
def clean(image):
    # normalized gamma curve
    gamma = sk.exposure.adjust_gamma(image, 1 / 2.2)
    gamma = renorm(gamma)
    # choose appropriate threshold for image
    gamma_filtered = ndi.median_filter(gamma, size=5)
    yen = sk.filters.threshold_yen(gamma_filtered)
    # mask background and small objects
    foreground = gamma_filtered >= yen
    # remove any objects less than or equal to ~0.1% area
    size = image.size // 1000
    foreground = sk.morphology.remove_small_objects(foreground, max_size = size)
    gamma[~foreground] = 0

    return gamma


# given an image and threshold, looks to segment out the galaxy
# does so by finding centermost object in the frame, them asking others
def find_galaxy(image, deblend=False, threshold=0):
    # make a segmentation map and catalog of objects in frame
    size = max(np.array(image.shape)) // 4
    cat, segm, segm_deblend = pf.make_catalog(
        image, plot=False, deblend=deblend, npixels=size, threshold=threshold
    )
    if deblend:
        segm = segm_deblend
    # the galaxy is the center most object in the frame save it
    if cat.nlabels == 1:
        galaxy = cat[0]
        center = True
    elif cat.nlabels == 0:
        print("Segmentation failed to find any regions.")
        return image, segm, False
    else:
        center_row, center_col = np.array(image.shape) // 2
        includes_center = [
            obj for obj in cat
            if obj.bbox.ixmin <= center_col <= obj.bbox.ixmax
            and obj.bbox.iymin <= center_row <= obj.bbox.iymax
        ]
        center = len(includes_center) > 0
        if center:
            galaxy = max(includes_center, key=lambda c: c.area)
        else:
            dist_from_center = np.array(
                list(
                    map(
                        lambda obj: np.linalg.norm(obj.centroid - [center_col, center_row]),
                        cat
                    )
                )
            )
            galaxy = min(cat, key=lambda obj: dist_from_center[obj.label - 1])
    return galaxy, segm, center


class Pipeline:
    def __init__(self, data_dir=Path.cwd(), seed=253, l=2**16, downscale=2):
        root = Path(data_dir)
        self.DATA = root / "galaxies_training"
        self.GRAY = root / "galaxies_grayscale"
        self.CROP = root / "galaxies_cropped"
        self.CLEAN = root / "galaxies_cleaned"
        self.RNG = np.random.default_rng(seed)
        self.SCALE = downscale
        self.L = l


    # Clean and rescale galaxies, saving gray and cropped images
    def process_images(self, limit=0):
        image_files = list(Path(self.DATA).glob("*.jpg"))
        Path(self.GRAY).mkdir(parents=True, exist_ok=True)
        Path(self.CROP).mkdir(parents=True, exist_ok=True)
        Path(self.CLEAN).mkdir(parents=True, exist_ok=True)
        # Process each file in the directory
        i = 0
        failed_images = list()
        for file in image_files:
            with warnings.catch_warnings():
                warnings.filterwarnings("error", category=RuntimeWarning)
                try:
                    # Load image, downscale to speed things along
                    save_name = file.stem+".png"
                    original = sk.io.imread(file, as_gray=True)
                    # Save the grayscale and cleaned image
                    cleaned = clean(original)
                    #sk.io.imsave(self.GRAY/save_name, sk.util.img_as_uint(cleaned))
                    # Save the cropped image of the galaxy
                    corrected, edge, center = self.correct(cleaned, original)
                    sk.io.imsave(self.CLEAN / save_name, sk.util.img_as_uint(corrected))
                    # add to failed list if either not centered or extended past border
                    if edge:
                        print(f"Galaxy {file.stem} crop extends to end of image")
                        failed_images.append(file.name)
                    if not center:
                        print(f"Could not find object in center for Galaxy {file.stem}")
                        failed_images.append(file.name)
                    i = i + 1
                    if i % 100 == 0:
                        print(f"Processed {i} images.")
                    if limit != 0 and i >= limit:
                        break
                except Exception as e:
                    print(f"Skipping {file.name}: {e}")
                    failed_images.append(file.name)
                    continue

        csv_path = Path(self.DATA) / "../failed_images.csv"
        with open(csv_path, "w") as f:
            f.write("filename\n")
            for name in failed_images:
                f.write(f"{name}\n")
        return failed_images


    # Clean and rescale galaxies, saving gray and cropped images
    def process(self, galaxy_id):
        Path(self.GRAY).mkdir(parents=True, exist_ok=True)
        Path(self.CROP).mkdir(parents=True, exist_ok=True)
        Path(self.CLEAN).mkdir(parents=True, exist_ok=True)
        file_name = galaxy_id+".jpg"
        save_name = galaxy_id+".png"
        original = sk.io.imread(self.DATA/file_name, as_gray=True)
        # Save the grayscale and cleaned image
        cleaned = clean(original)
        sk.io.imsave(self.GRAY/save_name, sk.util.img_as_uint(cleaned))
        # Save the cropped image of the galaxy
        corrected, edge, center = self.correct(cleaned, original)
        sk.io.imsave(self.CLEAN / save_name, sk.util.img_as_uint(corrected))


    # Finds the galaxy in the image, masks others, and returns cropped image
    def correct(self, image, original):
        # find the galaxy
        corrected = image.copy()
        galaxy, segm, center = find_galaxy(corrected)
        # mask everything that is not the galaxy
        mask = np.zeros(image.shape, dtype=bool)
        mask[segm.data == galaxy.label] = True
        corrected[~mask] = 0
        # correct for the background sky
        corrected = self.fill_sky(corrected, original[segm.data == 0])
        # crop down to just the galaxy itself
        bbox = galaxy.bbox
        padding = (max(np.array(image.shape)) // 40)
        margin = (
            max(
                min(bbox.ixmin, bbox.iymin, image.shape[1] - bbox.ixmax, image.shape[0] - bbox.iymax) - padding,
                0
            )
        )
        cropped = sk.util.crop(corrected, margin)
        size = max(np.array(original.shape)) // self.SCALE
        corrected = sk.transform.resize(cropped, (size, size),
                                        order=2, anti_aliasing=True)
        return corrected, cropped.shape == image.shape, center


    # adds back average sky based on statistics of original image, but smoothed at boundary
    # takes cleaned 2d array of image, and original, and pixel width of the taper
    def fill_sky(self, image, sky_pixels):
        # find original sky statistics
        mu, sigma = float(np.mean(sky_pixels)), float(np.std(sky_pixels))
        # Soft mask: Distance transform taper, 1%*size pixels wide
        taper = (max(np.array(image.shape)) // 100) + 1
        galaxy_mask = np.zeros_like(image, dtype=float)
        galaxy_mask[image>0] = 1.0
        soft_mask = ndi.gaussian_filter(galaxy_mask, sigma=taper / 2.0)
        soft_mask = renorm(soft_mask)
        # Synthetic sky, zero it where galaxy is
        sky = self.RNG.normal(mu, sigma, size=image.shape)
        # apply to galaxy, inside the galaxy remains untouched since sky is 0
        sky = (1.0 - soft_mask) * sky
        corrected = (image * soft_mask) + sky
        return corrected


    # Given a galaxy ID, plot the original, the cleaned, and the standardized versions
    def plot(self, galaxy_id):
        file_name = galaxy_id+".jpg"
        galaxy = sk.io.imread(self.DATA/file_name, as_gray=True)
        cleaned = clean(galaxy)
        corrected, edge, center = self.correct(cleaned, galaxy)
        fig, ax = plt.subplots(1, 5, figsize=(15, 3))
        ax[0].imshow(sk.io.imread(self.DATA/file_name))
        ax[0].axis("off")
        ax[1].imshow(sk.exposure.adjust_gamma(galaxy, 1 / 2.2), cmap="gray")
        ax[1].axis("off")
        ax[2].imshow(cleaned, cmap="gray")
        ax[2].axis("off")
        ax[3].imshow(corrected, cmap="gray")
        ax[3].axis("off")
        F, _, _= fft_pair(corrected)
        ax[4].imshow(log_scale(np.fft.fftshift(F)), cmap='gray', vmin=0, vmax=1)
        ax[4].axis("off")
        plt.show()
