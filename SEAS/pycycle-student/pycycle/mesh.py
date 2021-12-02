import numpy as np


class LineElement:
    """LineElement represents a finite arc in the discretization of the domain boundary."""
    def __init__(self, a, b, n, is_fault):
        """Constructor.

        :param a: Start point
        :param b: End point
        :param n: Outward-pointing normal
        :param is_fault: Flag this element as fault
        """
        self.a = np.array(a)
        self.h = np.array(b) - self.a
        self.h_norm = np.linalg.norm(self.h)
        self.n = np.array(n)
        self.n = self.n / np.linalg.norm(self.n)
        self.is_fault = is_fault

    def xi(self, theta):
        """Map from interval [-1, 1] to line a-b.

        :param theta: Scalar in [-1, 1].
        """
        # TODO: implement - DONE
        #return np.array([0.0, 0.0])
        return self.h * (theta + 1) / 2 + self.a

    def basis(self, theta):
        """Basis function evaluated at theta.

        :param theta: Scalar in [-1, 1]
        """
        # TODO: implement (constant in the finite case) - DONE
        #return 0.0
        return 1

    def factor(self, theta):
        """Integration factor.
        Must return basis(theta) * |xi'(theta)|

        :param theta: Scalar in [-1, 1]
        """
        # TODO: implement - DONE
        #return 0.0
        #|xi'(theta)| = |self.h / 2| = self.h_norm / 2
        return self.h_norm / 2
        

    def collocation_point(self):
        """Returns midpoint of line."""
        return self.xi(0)

    def __repr__(self):
        return 'LineElement({}, {})'.format(self.a, self.a + self.h)


class InfiniteLineElement:
    """InfiniteLineElement represents an infinite arc in the discretization of the domain boundary."""
    def __init__(self, a, n):
        """Constructor.

        :param a: Start point (line direction is also a)
        :param n: Outward-pointing normal
        :param is_fault: Flag this element as fault
        """
        self.a = np.array(a)
        self.a_norm = np.linalg.norm(self.a)
        self.n = np.array(n)
        self.n = self.n / np.linalg.norm(self.n)
        self.is_fault = False

    def xi(self, theta):
        """Map from interval [-1, 1] to line starting at "a" with direction "a" extending to infinity.

        :param theta: Scalar in [-1, 1].
        """
        # TODO: implement - DONE
        #return np.array([0.0, 0.0])
        
        return self.a*(theta+3)/(1-theta)
        # Why does the following not work?
        #constant p(theta)
        #return self.a * (-2)/(theta - 1)
        # or with steeper slope
        #return self.a * (4 + 2*theta)/(1 - theta)

    def basis(self, theta):
        """Basis function evaluated at theta.

        :param theta: Scalar in [-1, 1]
        """
        # TODO: implement - DONE
        # The basis function is $\frac{|a|^2}{|\xi|^2}$
        #return 0.0 
        #bas = self.a_norm**2 / (self.a_norm *(theta+3)/(1-theta))**2  --> simplified to
        bas = (1-theta)**2 / (theta+3)**2
        return bas


    def factor(self, theta):
        """Integration factor.
        Must return basis(theta) * |xi'(theta)|

        :param theta: Scalar in [-1, 1]
        """
        # TODO: implement - DONE
        #return 0.0
        bas = (1-theta)**2 / (theta+3)**2
        return bas * self.a_norm * 4 / (theta-1)**2
    

    def collocation_point(self):
        """Returns start point of line."""
        return self.xi(-1)

    def __repr__(self):
        return 'InfiniteLineElement({})'.format(self.a)


def tessellate_line(a, b, resolution, normal, is_fault=False):
    """Tessellate the line from a to b into small arcs, such that
       the arc length is smaller than resolution.

    :param a: Start point
    :param b: End point
    :param resolution: Target arc length
    :param normal: Outward-pointing normal
    :param is_fault: Flag all line elements as fault
    """
    origin = np.array(a)
    h = np.array(b) - origin
    N = int(np.ceil(np.linalg.norm(h) / resolution))
    return [
        LineElement(origin + n / N * h, origin + (n + 1) / N * h, normal,
                    is_fault) for n in range(N)
    ]


def num_fault_elements(mesh):
    """Counts number of fault elements in mesh.

    :param mesh: List of LineElements.
    """
    count = 0
    for m in mesh:
        if m.is_fault:
            count = count + 1
    return count


def line_normal(a, b, star_centre):
    """Computes the outward-pointing of a line from a to b.
       The parameter star_centre here defines the inside,
       i.e. dot(star_centre-a, normal) < 0.

    :param a: Start point
    :param b: End point
    :param star_centre: Star centre of a star-shaped domain
    """
    c = np.array(star_centre) - np.array(a)
    normal = np.array(b) - np.array(a)
    normal[0], normal[1] = -normal[1], normal[0]
    normal /= np.linalg.norm(normal)
    return normal if np.inner(c, normal) < 0 else -normal
