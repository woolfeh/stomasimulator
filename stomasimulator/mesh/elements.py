""" Classes for 8-node brick and 6-node prismatic elements """

import numpy as np

from stomasimulator.geom.point import Point


class BlockElement(object):
    """ Base class for the 3D elements """

    element_id = 1

    def __init__(self, nodes):
        """
        It is assumed that the list is ordered so that (1-based) indices...

        hex8:
            1-4-3-2 and 5-6-7-8 form two separate planes with (i+4)th node above the ith node (e.g. 5 above 1)
        penta6:
            1-3-2 and 4-5-6 for two separate plane with (i+3)th node above the ith node (e.g. 4 above 1)

        :param nodes: list of Point objects or integer point ids
        :type nodes: list
        :return:
        """

        if isinstance(nodes, (list, tuple)) and len(nodes) == self.num_nodes():
            if isinstance(nodes[0], Point):
                self.node_ids = [n.id for n in nodes]
            else:
                self.node_ids = list(nodes)

            self.fibre_vector = None
            self.reference_element = None

            self.id = BlockElement.element_id
            BlockElement.element_id += 1
        else:
            raise ValueError('Error: must pass a list of length {} to create a {} object'
                             .format(self.num_nodes(), self.element_type()))

    def __repr__(self):
        return '{:4} - {:8}: {}'.format(self.id, self.element_type(), self.node_ids)

    @property
    def formatted_node_list(self):
        """ Get comma separated list of *global* node ids """
        return ''

    def facets(self):
        """ Get list of facets defined using the *global* node ids """
        pass

    def get_edges(self):
        """ Get list of edges defined using the *global* node ids """
        edges = []

        for le in self.local_edges():
            edges.append((self.node_ids[le[0] - 1], self.node_ids[le[1] - 1]))

        return edges

    @classmethod
    def element_type(cls):
        """ Get the elements type - returns a string """
        pass

    @classmethod
    def num_nodes(cls):
        """ Get the number of nodes for the element """
        pass

    @classmethod
    def local_facets(cls):
        """ Get list of facets defined using local node numbers """
        pass

    @classmethod
    def local_edges(cls):
        """ Get list of edges using local node numbers """
        return []

    @classmethod
    def build_element_from_node_ids(cls, node_ids):
        """
        Factory method that will construct an element from a list of node ids
        :param node_ids: list of integer node ids
        :return: Penta6Element or Hex8Element
        """
        if len(node_ids) == Penta6Element.num_nodes():
            return Penta6Element(node_ids)
        elif len(node_ids) == Hex8Element.num_nodes():
            return Hex8Element(node_ids)
        else:
            raise NotImplementedError("Shouldn't be here! an element with {} node(s) is unsupported"
                                      .format(len(node_ids)))

    def jacobian_det(self, nodes_map, r=0.0, s=0.0, t=0.0):
        """
        Get the Jacobian determinant
        :param nodes_map:
        :param r:
        :param s:
        :param t:
        :return:
        """
        pass

    def jacobian_ratio(self, nodes_map):
        """
        Calculate the ratio of the max nodal Jacobian to the min nodal Jacobian.
        itertools is used to generate the parametric nodal coordinates [ [-1,-1,-1], [-1,-1,1], etc. ]
        ANSYS: '...a Jacobian ratio >40 computed at element nodal points will cause a mesh failure"
        """

        from itertools import product
        jacs = [self.jacobian_det(nodes_map, r=a[0], s=a[1], t=a[2]) for a in product([-1.0, 1.0], repeat=3)]

        return max(jacs) / min(jacs)


####################################################################################################

class Hex8Element(BlockElement):
    """ Represents an 8-node brick element """

    @property
    def formatted_node_list(self):
        return '{},{},{},{},{},{},{},{}'.format(*self.node_ids)

    @classmethod
    def element_type(cls):
        return 'hex8'

    @classmethod
    def num_nodes(cls):
        return 8

    @classmethod
    def local_facets(cls):
        return (1, 2, 6, 5), \
               (2, 3, 7, 6), \
               (3, 4, 8, 7), \
               (1, 5, 8, 4), \
               (1, 4, 3, 2), \
               (5, 6, 7, 8)

    def facets(self):
        return [(self.node_ids[lf[0] - 1],
                 self.node_ids[lf[1] - 1],
                 self.node_ids[lf[2] - 1],
                 self.node_ids[lf[3] - 1])
                for lf in Hex8Element.local_facets()]

    @classmethod
    def local_edges(cls):
        return (1, 2), \
               (1, 4), \
               (1, 5), \
               (2, 3), \
               (2, 6), \
               (3, 4), \
               (3, 7), \
               (4, 8), \
               (5, 6), \
               (5, 8), \
               (6, 7), \
               (7, 8)

    def jacobian_det(self, nodes_map, r=0.0, s=0.0, t=0.0):
        rm, rp = 1 - r, 1 + r
        sm, sp = 1 - s, 1 + s
        tm, tp = 1 - t, 1 + t

        shape = [[-sm * tm, sm * tm, sp * tm, -sp * tm, -sm * tp, sm * tp, sp * tp, -sp * tp],
                 [-rm * tm, -rp * tm, rp * tm, rm * tm, -rm * tp, -rp * tp, rp * tp, rm * tp],
                 [-rm * sm, -rp * sm, -rp * sp, -rm * sp, rm * sm, rp * sm, rp * sp, rm * sp]]

        pts = [nodes_map[nid].xyz for nid in self.node_ids]

        # get matrix product - divide by 8 now because the shape matrix should have been divided by 8
        jacobian = np.dot(shape, pts) / 8

        return np.linalg.det(jacobian)


class Penta6Element(BlockElement):
    """ Represents a 6-node pentahedral element """

    @property
    def formatted_node_list(self):
        return '{},{},{},{},{},{}'.format(*self.node_ids)

    @classmethod
    def element_type(cls):
        return 'penta6'

    @classmethod
    def num_nodes(cls):
        return 6

    @classmethod
    def local_facets(cls):
        return (1, 3, 2), (1, 2, 5, 4), (2, 3, 6, 5), (1, 4, 6, 3), (4, 5, 6)

    def facets(self):
        """ Get list of facets """

        facet_list = []
        for lf in Penta6Element.local_facets():
            facet_list.append(self.global_facet(lf))

        return facet_list

    def global_facet(self, local_facet):
        """
        Get list of global node ids that represent the facet
        :param local_facet:
        :return:
        """
        return [self.node_ids[local_nid - 1] for local_nid in local_facet]

    @classmethod
    def local_edges(cls):
        return (1, 2), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6), (4, 5), (4, 6), (5, 6)

    def jacobian_det(self, nodes_map, r=0.0, s=0.0, t=0.0):
        shape = [[-(1 - t), 1 - t, 0.0, -(1 + t), 1 + t, 0.0],
                 [-(1 - t), 0.0, 1 - t, -(1 + t), 0.0, 1 + t],
                 [-(1 - r - s), -r, -s, 1 - r - s, r, s]]

        pts = [nodes_map[nid].xyz for nid in self.node_ids]

        # get matrix product - divide by 2 now because the shape matrix should have been divided by 2
        jacobian = np.dot(shape, pts) / 2

        return np.linalg.det(jacobian)


if __name__ == '__main__':
    pass
