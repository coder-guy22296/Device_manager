#!/usr/bin/python
import sys
import os
sys.path.append(os.getcwd())

import inLib


class Control(inLib.Module):

    def __init__(self, control, settings):
        print 'Initializing HotSpots.'
        inLib.Module.__init__(self, control, settings)
        self._hotspots = []


    def getHotSpots(self, optimize=False, priority=None):
        ''' A list of the current hot spots. '''
        if optimize:
            return self._optimize_order(priority)
        else:
            return self._hotspots


    def add(self):
        '''
        Adds a new hotspot. The HotSpot module will record the position of
        all relevant devices, as listed in the settings file.
        '''
        coordinates = {}
        for device in self._settings:
            control = inLib.get_nested_attr(self._control, device)
            get_position = inLib.get_nested_attr(control,
                                                 self._settings[device][0])
            coordinates[device] = get_position()
        self._hotspots.append(HotSpot(coordinates, self._settings, self._control))
        return self._hotspots[-1]


    def remove(self, hotspot):
        '''
        Removes a previously recorded hot spot from the :attr:`hotSpots` list.

        :Parameters:
            *hotspot*:
                A :class:`HotSpot` instance.
        '''
        self._hotspots.remove(hotspot)


    def _optimize_order(self, device):
        # A list with optimized order:
        optimizedList = []
        remaining = [_ for _ in self._hotspots]
        # We start with the last one:
        optimizedList.append(remaining.pop())

        # Figure out the number of coordinates:
        coordinates = optimizedList[-1].getCoordinates()[device]
        try:
            nCoordinates = len(coordinates)
        except TypeError:
            nCoordinates = 1
        except:
            raise

        while len(remaining) > 0:
            # A list of distances of all remaining hotspots to the last of the
            # optimized list:
            distances = []
            coordinates_last = optimizedList[-1].getCoordinates()[device]
            # Loop through remaining hotspots to fill up distances list:
            for hotspot in remaining:
                coordinates_this = hotspot.getCoordinates()[device]
                if nCoordinates > 1:
                    # If there are multiple coordinates, we want to choose the one
                    # with maximal distance.
                    ds = []
                    for _ in range(nCoordinates):
                        d = abs(coordinates_last[_] - coordinates_this[_])
                        ds.append(d)
                    d = max(ds)
                else: 
                    d = abs(coordinates_last - coordinates_this)
                distances.append(d)
            # The minimal distance to the next hotspots:
            min_d = min(distances)
            # The index of the next hotspot:
            min_index = distances.index(min_d)
            # Append the next hotspot the optimized list:
            optimizedList.append(remaining.pop(min_index))
        return optimizedList
                
            



class HotSpot:
    def __init__(self, coordinates, settings, control):
        self._coordinates = coordinates
        self._settings = settings
        self._control = control

    def getCoordinates(self):
        return self._coordinates

    def moveTo(self):
        '''
        Moves to the  hot spot. Specifically, every relevant device, as listed
        in the settings file, will move to the previously recorded position.
        '''
        for device in self._settings:
            control = inLib.get_nested_attr(self._control, device)
            set_position = inLib.get_nested_attr(control,
                                                 self._settings[device][1])
            set_position(self._coordinates[device])


if __name__ == '__main__':
    control = Control(None, {})
    control._hotspots = [
            {'hal4000' : [200,300,0], 'scope.PFS' : 200},
            {'hal4000' : [100,200,0], 'scope.PFS' : 100},
            {'hal4000' : [133,421,0], 'scope.PFS' : 800}
            ]
    print control.getHotSpots()
    control.optimizeOrder('hal4000')
    print control.getHotSpots()
