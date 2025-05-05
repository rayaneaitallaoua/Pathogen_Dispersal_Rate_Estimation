import os

def set_tests_dir():

    # Check if we are already in 'Tests' directory
    if os.path.basename(os.getcwd()) == "Tests":
        print("Already in the 'Tests' directory!")
        print(f"Current working directory: {os.getcwd()}")
    else:
        # Check if 'Tests' exists in current directory
        if not os.path.exists("Tests"):
            os.makedirs("Tests")
        os.chdir("Tests")
        print(f"Moved into 'Tests' directory. Current working directory: {os.getcwd()}")


import matplotlib.pyplot as plt


import matplotlib.pyplot as plt

def plot_sampling(lattice_size_x,
                               lattice_size_y,
                               sampled_positions,
                               r,
                               title='Sampling Visualization'):
    """
    Visualize a lattice grid and highlight sampled nodes with a sampling circle.

    Parameters:
    - lattice_size_x: width of the grid
    - lattice_size_y: height of the grid
    - sampled_positions: list or set of (x, y) tuples representing sampled nodes
    - r: sampling radius (used to draw the reference circle)
    - title: title of the plot
    """
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_aspect('equal')
    ax.set_xlim(-1, lattice_size_x)
    ax.set_ylim(-1, lattice_size_y)

    # Lattice center
    center_x, center_y = lattice_size_x // 2, lattice_size_y // 2

    # Draw all lattice positions
    for x in range(lattice_size_x):
        for y in range(lattice_size_y):
            color = 'lightgray'
            edge = 'black'
            lw = 0.5
            if (x, y) in sampled_positions:
                color = 'red'
                lw = 1.0
            circle = plt.Circle((x, y), 0.3, color=color, edgecolor=edge, lw=lw)
            ax.add_patch(circle)

    # Mark center
    ax.plot(center_x, center_y, 'ko', markersize=6, label='Center')

    # Draw sampling circle
    sampling_circle = plt.Circle((center_x, center_y), r, color='blue', fill=False, linestyle='--', linewidth=1.5, label=f'r = {r}')
    ax.add_patch(sampling_circle)

    # Decorations
    ax.set_xticks(range(lattice_size_x))
    ax.set_yticks(range(lattice_size_y))
    ax.grid(True, linestyle='--', alpha=0.3)
    ax.set_title(title)
    plt.xlabel("X")
    plt.ylabel("Y")
    ax.legend(loc='upper right')
    plt.tight_layout()
    plt.show()