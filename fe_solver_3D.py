import numpy as np
# import matplotlib.pyplot as plt

nodes_test = [(0, 0, 0, 0, 0, 0), (5, 0, 0, 1, 1, 0), (0, 5, 0, 0, 0, 0)] #0 = restrained
forces_test = [(1, 0, -5, 0)] #node number, x_value, y_value
E_test = 200e6 #kN/m2
A_test = 100e-7
bars_test = [(0, 1), (1, 2)]

class FESolver():
    def __init__(self, nodes, forces, bars, E, A):
        self.nodes = nodes
        self.forces = forces
        self.E = E
        self.A = A
        self.fe_bars = bars


    class Bar():
        def __init__(self, outer_instance, n_1, n_2, A):
            self.outer_instance = outer_instance
            self.n_1 = n_1
            self.n_2 = n_2
            self.node1 = self.outer_instance.nodes[n_1][:3]
            self.node2 = self.outer_instance.nodes[n_2][:3]
            self.A = A
            
            self.vec = [self.node2[0] - self.node1[0], self.node2[1] - self.node1[1], self.node2[2] - self.node1[2]]
            self.length = np.sqrt(self.vec[0]**2 + self.vec[1]**2 + self.vec[2]**2)
            self.cosines = [self.vec[0]/self.length, self.vec[1]/self.length, self.vec[2]/self.length]
            
            self.T_m = np.matrix([[self.cosines[0], 0], 
                                    [self.cosines[1], 0],
                                    [self.cosines[2], 0],
                                    [0, self.cosines[0]],
                                    [0, self.cosines[1]],
                                    [0, self.cosines[2]]
                                    ]).round(4)
            
            K_m_local = np.matrix([[1, -1], [-1, 1]])
            self.K = self.T_m*(self.A*self.outer_instance.E/self.length)*K_m_local*(self.T_m.T)
        
    def solve(self):
        self.bars = []
        for bar in self.fe_bars:
            self.bars.append(self.Bar(self, bar[0], bar[1], self.A))

        matrix_size = 3*len(self.nodes)

        delta = np.zeros((matrix_size, ))
        for i, node in enumerate(self.nodes):
            delta[3*i] = node[3]
            delta[3*i + 1] = node[4]
            delta[3*i + 2] = node[5]

        q = np.zeros((matrix_size, ))
        for i, f in enumerate(self.forces):
            q[3*f[0]] = f[1]
            q[3*f[0] + 1] = f[2]
            q[3*f[0] + 2] = f[3]

        global_matrix = np.matrix(np.zeros((matrix_size, matrix_size)))

        for bar in self.bars:
            i, j = bar.n_1, bar.n_2
            global_matrix[3*i: 3*i+3, 3*i: 3*i+3] = global_matrix[3*i: 3*i+3, 3*i: 3*i+3] + bar.K[0:3, 0:3]
            global_matrix[3*j: 3*j+3, 3*j: 3*j+3] = global_matrix[3*j: 3*j+3, 3*j: 3*j+3] + bar.K[3:6, 3:6]
            global_matrix[3*j: 3*j+3, 3*i: 3*i+3] = global_matrix[3*j: 3*j+3, 3*i: 3*i+3] + bar.K[3:6, 0:3]
            global_matrix[3*i: 3*i+3, 3*j: 3*j+3] = global_matrix[3*i: 3*i+3, 3*j: 3*j+3] + bar.K[3:6, 0:3]

        global_matrix = global_matrix.round(2)

        free_dofs = []
        specified_dofs = []
        for i, dis in enumerate(delta):
            if dis == 1:
                free_dofs.append(i)
            else:
                specified_dofs.append(i)

        rearranged_matrix = global_matrix.copy()
        rearranged_matrix[0:len(free_dofs), 0:len(free_dofs)] = global_matrix[free_dofs, :][:, free_dofs]
        rearranged_matrix[len(free_dofs):matrix_size+1, len(free_dofs):matrix_size+1] = global_matrix[specified_dofs, :][:, specified_dofs]
        rearranged_matrix[len(free_dofs):matrix_size+1, 0:len(free_dofs)] = global_matrix[specified_dofs, :][:, free_dofs]
        rearranged_matrix[0:len(free_dofs), len(free_dofs):matrix_size+1] = global_matrix[free_dofs, :][:, specified_dofs]

        u = np.linalg.solve(rearranged_matrix[0:len(free_dofs), 0:len(free_dofs)], q[free_dofs])

        displacements = np.zeros(matrix_size)
        displacements[specified_dofs] = delta[specified_dofs]
        displacements[free_dofs] = u

        reactions = global_matrix * np.matrix(displacements.reshape(-1, 1)).round(3)
            
        ratio_x = 0.15*(max([node[0] for node in self.nodes]) - min([node[0] for node in self.nodes]))/max([abs(x) for x in displacements[0::2]])
        ratio_y = 0.15*(max([node[1] for node in self.nodes]) - min([node[1] for node in self.nodes]))/max([abs(y) for y in displacements[1::2]])
         
        ratio = min([ratio_x, ratio_y])

        self.plotting_nodes = []
        for i, node in enumerate(self.nodes):
            self.plotting_nodes.append((node[0] + ratio*displacements[3*i], node[1] + ratio*displacements[3*i+1], node[2] + ratio*displacements[3*i+2]))

        

solver = FESolver(nodes_test, forces_test, bars_test, E_test, A_test)
solver.solve()

# for bar in bars:
#     plt.plot([bar.node1[0], bar.node2[0]], [bar.node1[1], bar.node2[1]], 'k-')
#     plot_displacements = 50*displacements
#     delta_x1, delta_y1 = bar.node1[0] + plot_displacements[2*(bar.n_1)], bar.node1[1] + plot_displacements[2*(bar.n_1)+1]
#     delta_x2, delta_y2 = bar.node2[0] + plot_displacements[2*(bar.n_2)], bar.node2[1] + plot_displacements[2*(bar.n_2)+1]
#     print((delta_x1, delta_y1), (delta_x2, delta_y2))
#     plt.plot([delta_x1, delta_x2], [delta_y1, delta_y2], '--', color='blue')