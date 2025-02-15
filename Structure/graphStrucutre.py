import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_binding_locations(df):
    # Create a 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Extract coordinates and binding locations
    x = df['Coordinates'].apply(lambda coord: coord[0])
    y = df['Coordinates'].apply(lambda coord: coord[1])
    z = df['Coordinates'].apply(lambda coord: coord[2])
    labels = df['Binding_Location']

    # Plot points in 3D
    ax.scatter(x, y, z, color='b', marker='o')

    # Annotate each point with its binding location
    for i, label in enumerate(labels):
        ax.text(x.iloc[i], y.iloc[i], z.iloc[i], str(label), size=8, zorder=1)

    # Set axis labels
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_zlabel('Z Coordinate')

    # Show plot
    plt.title(f"Binding Locations for {df['KOID'][0]}")
    plt.show()