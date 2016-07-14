#!/usr/bin/python


import inLib


class UI(inLib.ModuleUI):
    
    def __init__(self, control, ui_control):
        print 'Initializing HotSpots UI.'
        inLib.ModuleUI.__init__(self, control, ui_control,
                                'modules.hotSpots.hotSpots_design')

        self._ui.pushButtonAdd.clicked.connect(self.add)
        self._ui.pushButtonRemove.clicked.connect(self.remove)

        self._ui.listWidgetHotSpots.itemActivated.connect(self._on_row_activated)


    def add(self):
        hotspot = self._control.add()
        self._ui.listWidgetHotSpots.addItem(str(hotspot.getCoordinates()))
        items = []
        for i in xrange(self._ui.listWidgetHotSpots.count()):
            items.append(self._ui.listWidgetHotSpots.item(i))
        self._ui.listWidgetHotSpots.setItemSelected(items[-1], True)


    def _on_row_activated(self, item):
        hotspot = self._item2hotspot(item)
        hotspot.moveTo()


    def remove(self):
        item = self._ui.listWidgetHotSpots.selectedItems()[0]
        hotspot = self._item2hotspot(item)
        row = self._ui.listWidgetHotSpots.row(item)
        self._ui.listWidgetHotSpots.takeItem(row)
        self._control.remove(hotspot)


    def _item2hotspot(self, item):
        coordinatestring = str(item.text())
        hotspots = self._control.getHotSpots()
        for hotspot in hotspots:
            if coordinatestring == str(hotspot.getCoordinates()):
                return hotspot

