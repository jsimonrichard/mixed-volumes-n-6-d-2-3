from math import prod


def property_memoize(func):
    cache = {}

    def get(self):
        key = id(self)
        if key not in cache:
            cache[key] = func(self)
        return cache[key]

    return property(get, None, None, doc=func.__doc__)


default_bodies = [
    Polyhedron([[0, 0], [1, 0]]),
    Polyhedron([[0, 0], [0, 1]]),
    Polyhedron([[0, 0], [1, 1]]),
]


def mixed_volume(A: Polyhedron, B: Polyhedron):
    """
    Computes the mixed volume (times 2) of the two polyhedra A and B.
    """
    result = (A + B).volume() - A.volume() - B.volume()
    # Cast to the base ring of A to avoid weird bug that changes
    # volume output to Real Double Field
    return A.base_ring()(result)


class MVolMapProblem(object):
    """
    This class saves all of the parameters related to one version
    of the mixed volume map problem.
    """

    def __init__(
        self,
        factors: int,
        mixed_volume=mixed_volume,
        bodies=default_bodies,
    ):

        self.dimension = 2 # these code probably doesn't work for other dimensions (unchecked)
        self.factors = factors
        self.n = self.dimension * factors
        self.mixed_volume = mixed_volume
        self.bodies = bodies

    @property_memoize
    def coord_partitions(self):
        """
        We index products of mixed volumes (the coordinates of the mvol_map
        output) by partitions. As one example, if we are in dimension
        d=2 and we have 3 factors in the product of mixed volumes, then
        we have {{0,1},{2,3},{4,5}} -> V(K0,K1)V(K2,K3)V(K4,K5) as
        one of the 15 partitions (products).

        d = dimension = number of elements in a part
        num_parts = number of parts

        number of bodies is d * num_parts

        To access this property, use `<problem>.coord_partitions`.
        """
        # See https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/set_partition.html
        return SetPartitions(self.n, [self.dimension] * self.factors)

    def mvol_map(self, body_config, mixed_volume=None):
        """
        Computes our mixed volume_ map for a list of bodies.
        """

        if mixed_volume is None:
            mixed_volume = self.mixed_volume

        return vector(
            [
                prod(
                    mixed_volume(*map(lambda i: body_config[i - 1], mv_group))
                    for mv_group in part
                )
                for part in self.coord_partitions
            ]
        )

    @property_memoize
    def configuration_partitions(self):
        """
        We index configurations of bodies by partitions.

        We are only interested in configurations for which all bodies have
        multiplicities greater than two (to remove linearity) and less than
        factors (to prevent all terms collapsing to zero); the partitions are
        limited accordingly.

        If the order of the bodies (before a configuration is applied) matters,
        set ordered=True.

        To access this property, use `<problem>.configuration_partitions`.
        """

        partitions = []
        for seg_comb_unordered in Partitions(self.n, min_part=2, max_part=self.factors):
            if len(seg_comb_unordered) <= 3:
                # Any 3 (or less) segments can be transformed into any other 3 (or less) segments
                partitions.extend(SetPartitions(self.n, seg_comb_unordered))
            else:
                for seg_comb in Permutations(seg_comb_unordered):
                    partitions.extend(
                        OrderedSetPartitions(self.n, Composition(seg_comb))
                    )

        return partitions

    def get_body_config(self, partition):
        """
        Constructs the configuration of bodies indexed by the partition.
        """
        configuration = [None] * partition.size()
        for i, part in enumerate(partition):
            for j in part:
                configuration[j - 1] = self.bodies[i]
        return configuration

    @property_memoize
    def image_rays(self):
        """
        Computes the extremes (rays) of mvol_map's image.

        To access this property, use `<problem>.image_rays`.
        """
        rays = []
        for partition in self.configuration_partitions:
            configuration = self.get_body_config(partition)
            rays.append(self.mvol_map(configuration))
        return rays

    @property_memoize
    def image_convex_hull(self):
        """
        Computes the convex hull of mvol_map's image.

        To access this property, use `<problem>.image_convex_hull`.
        """
        rays = self.image_rays
        return Polyhedron(rays=rays)

    @property_memoize
    def image_ineqs(self):
        """
        The inequalities of mvol_map's image.

        To access this property, use `<problem>.image_ineqs`.
        """
        return [vector(ineq[1:]) for ineq in self.image_convex_hull.inequalities_list()]

    @property_memoize
    def _coord_part_inv_map(self):
        """
        Computes the inverse of the coordinate partition map.

        To access this property, use `<problem>._coord_part_inv_map`.
        """
        return dict([[part, i] for i, part in enumerate(self.coord_partitions)])

    def permute_vec(self, vec: vector, perm: Permutation):
        """
        Applies the action of a permutation of bodies (in S_{d*factors}) on a
        vector in the image of mvol_map.
        """
        result = zero_vector(vec.base_ring(), len(self.coord_partitions))
        for i, x in enumerate(vec):
            part = permute_partition(self.coord_partitions[i], perm)
            j = self._coord_part_inv_map[part]
            result[j] = x
        return result

    def break_into_orbits(self, vectors):
        """
        Breaks a complete set (or iterator) of vectors in the image of mvol_map
        into orbits under the action of S_{d*factors}.

        This works because we already know all of the vectors (and therefore
        no additional vectors need to be added to the graph that is used to find
        the orbits).
        """

        G = Graph()
        for vec in vectors:
            # Allow hashing
            vec = vector(vec)  # copy to prevent side effects
            vec.set_immutable()

            G.add_vertex(vec)
            vec2 = self.permute_vec(vec, Permutation([2, 3, 4, 5, 6, 1]))
            vec2.set_immutable()
            vec3 = self.permute_vec(vec, Permutation([2, 1, 3, 4, 5, 6]))
            vec3.set_immutable()
            if vec != vec2:
                G.add_edge([vec, vec2])
            if vec != vec3:
                G.add_edge([vec, vec3])
        return G.connected_components()

    @property_memoize
    def image_ineq_orbits(self):
        """
        The orbits of the inequalities of mvol_map's image.

        To access this property, use `<problem>.image_ineq_orbits`.
        """
        return self.break_into_orbits(self.image_ineqs)

    @property_memoize
    def _mvol_map_config_part_inv_dict(self):
        configs_ray_dict = {}
        for part in self.configuration_partitions:
            configuration = self.get_body_config(part)
            ray = self.mvol_map(configuration)
            ray /= max(ray)  # Scale largest entry to 1
            configs_ray_dict[str(ray)] = part
        return configs_ray_dict

    def get_config_from_vec(self, vec):
        """
        Computes the inverse of the function
        configuration partition -> configuration -> mvol_map.
        """
        vec /= max(vec)  # Scale largest entry to 1
        return self._mvol_map_config_part_inv_dict[str(vec)]


def permute_partition(partition, perm: Permutation):
    """
    Applies a permutation of bodies (in S_{d*factors}) to a coordinate partition.
    """
    return SetPartition([[perm(x) for x in part] for part in partition])


def get_config_type(coord_part, ordered=False):
    """
    Returns the configuration type of a coordinate partition.
    """
    result = [len(s) for s in coord_part]
    # not sure if this is needed; only for regularizing sorted partitions (compositions) into unordered ones
    if ordered:
        result = sorted(result, reverse=True)
    return result
