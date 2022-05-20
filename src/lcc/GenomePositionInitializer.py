# from mpl_toolkits import mplot3d

# %matplotlib inline
import os, sys
import numpy as np
import matplotlib.pyplot as plt

'''
For initializing Genome location in a 3D cell container
'''


def DetermineDistance_3D(P1, P2):
    return np.sqrt(np.sum((P1 - P2) ** 2, axis=1))


def TrimToShape(Nodes, Shape):
    if Shape == 'rectangle':
        return Nodes
    else:
        print('Unrecognizable Shape for Trimming: ', Shape)
        print('Available Trimming Shapes: rectangle')
        sys.exit(0)


def GetNodesAndDistances(Dim_X, Dim_Y, Dim_Z, N_Nodes, shape='rectangle'):
    assert (isinstance(N_Nodes, int)) & (
                N_Nodes > 1), 'ERROR: N_Nodes (Input: %s) must be an integer greater than 1.' % N_Nodes

    XYZ = np.array([Dim_X, Dim_Y, Dim_Z])

    # Generate sets of three random numbers
    Nodes_Original = np.random.rand(N_Nodes, 3) * XYZ

    # TODO:Potential node trimming step for shape control
    Nodes_Trimmed = TrimToShape(Nodes_Original, shape)

    # Arrange Nodes by finding the nearest neighbor
    Node_Reference = Nodes_Original[0]
    Nodes_Modified = Nodes_Original[1:]  # Nodes excluding the reference node

    Nodes_Arranged = np.zeros_like(Nodes_Original)
    Nodes_Arranged[0] = Node_Reference

    Distances = np.zeros_like(Nodes_Original)

    # print('Prev', Nodes_Modified)
    # print('Node:', Node_Reference)

    for i in range(Nodes_Modified.shape[0]):
        Distance = DetermineDistance_3D(Node_Reference, Nodes_Modified)
        Idx_NearestNeighbor = Distance.argmin()
        Distances[i] = Distance[Idx_NearestNeighbor]
        Node_Reference = Nodes_Modified[Idx_NearestNeighbor]
        Nodes_Arranged[i + 1] = Node_Reference
        Nodes_Modified = np.delete(Nodes_Modified, Idx_NearestNeighbor, axis=0)
        # print(i)
        # print('Node    :', Node_Reference)
        # print('Rest    :', Nodes_Modified)
        # print('Arranged:', Nodes_Arranged)

    return Nodes_Arranged, Distances


def Plot3D(Nodes):
    ax = plt.axes(projection='3d')

    # Data for a three-dimensional line
    X = Nodes[:, 0]
    Y = Nodes[:, 1]
    Z = Nodes[:, 2]
    ax.plot3D(X, Y, Z, 'red')

    # Data for three-dimensional scattered points
    # ax.scatter3D(X, Y, Z, c='blue')
    # ax.scatter3D(X, Y, Z, c=Z, cmap='Greens')

    plt.show()


def main():   # add verbose
    Dim_X = 10000
    Dim_Y = 20000
    Dim_Z = 10000
    N_Nodes = 10000
    Nodes, Distances = GetNodesAndDistances(Dim_X, Dim_Y, Dim_Z, N_Nodes)
    Plot3D(Nodes)


if __name__ == '__main__':
    main()



