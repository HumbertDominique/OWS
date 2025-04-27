import numpy as np
import matplotlib.pyplot as plt

def simulate_lens(phase_screen, focal_length, wavelength, pixel_size):
    """
    Simulates the effect of a lens on an incoming wavefront.

    Parameters:
    - phase_screen: 2D numpy array representing the phase screen in the pupil plane (in radians)
    - focal_length: focal length of the lens (in meters)
    - wavelength: wavelength of the light (in meters)
    - pixel_size: physical size of each pixel in the pupil plane (in meters)

    Returns:
    - intensity: 2D array representing the intensity in the focal plane
    """

    # Grid size
    N = phase_screen.shape[0]
    k = 2 * np.pi / wavelength

    # Coordinate grids
    x = np.linspace(-N/2, N/2 - 1, N) * pixel_size
    X, Y = np.meshgrid(x, x)

    # Apply the lens: quadratic phase function
    lens_phase = np.exp(-1j * k * (X**2 + Y**2) / (2 * focal_length))

    # Input field (wavefront): phase screen + lens phase
    input_field = np.exp(1j * phase_screen) * lens_phase

    # Compute field at the focal plane using FFT (Fraunhofer approximation)
    focal_field = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(input_field)))

    # Intensity pattern (simulated image)
    intensity = np.abs(focal_field)**2

    return intensity

# Example usage
if __name__ == "__main__":
    # Simulation parameters
    N = 512  # grid size
    wavelength = 500e-9  # 500 nm
    focal_length = 0.2  # 20 cm
    pixel_size = 10e-6  # 10 microns

    # Generate a random phase screen for demonstration
    phase_screen = np.random.uniform(-np.pi, np.pi, (N, N))

    # Simulate
    intensity = simulate_lens(phase_screen, focal_length, wavelength, pixel_size)

    # Display result
    plt.imshow(intensity, cmap='hot', extent=[-1, 1, -1, 1])
    plt.title('Intensity at the Focal Plane')
    plt.colorbar(label='Intensity')
    plt.xlabel('x (arbitrary units)')
    plt.ylabel('y (arbitrary units)')
    plt.show()
