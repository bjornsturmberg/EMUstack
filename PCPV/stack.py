import layer

class Stack(object):
    """ Represents a stack of layers evaluated at one frequency.

        This includes the semi-infinite input and output layers.

        INPUTS:

          - `layers` : a tuple of :ThinFilm:s and :NanoStruct:s ordered
            from top to bottom layer.
    """
    def __init__(self, s_layers):
        self.layers = tuple(s_layers)

    def calc_scat(self):
        """ Calculate the transmission and reflection matrices of the stack"""
        for top, bottom in zip(self.layers[:-2], self.layers[1:]):
            pass

def god_function(raw_layers, lights, simmoargs = {}):
    for lay in layers:
        # Make the mesh (if appropriate)
        lay.make_mesh()

    stack_list = []
    for li in lights:
        # Find modes of the layer
        s_layers = (lay.eval(li, simmoargs) for lay in layers)
        stack = Stack(s_layers)
        stack.calc_scat()
        stack.delete_working()
