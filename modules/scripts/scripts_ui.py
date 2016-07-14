#!/usr/bin/python


from PyQt4 import QtCore, QtGui, Qsci
import inLib
import os



class UI(inLib.ModuleUI):
    
    def __init__(self, control, ui_control):
        print 'Initializing Scripts UI.'
        inLib.ModuleUI.__init__(self, control, ui_control, 'modules.scripts.scripts_design')

        self._runner = None

        self._ui.listWidgetFiles.dragEnterEvent = self._drag_enter_event
        self._ui.listWidgetFiles.dragMoveEvent = self._drag_move_event
        self._ui.listWidgetFiles.dropEvent = self._drop_event

        self._ui.listWidgetFiles.itemClicked.connect(self._on_file_clicked)
        self._ui.listWidgetFiles.itemActivated.connect(self.run)

        self._ui.pushButtonRun.clicked.connect(self.run)
        self._ui.pushButtonSave.clicked.connect(self._on_save_clicked)
        self._ui.pushButtonNew.clicked.connect(self.new)
        self._ui.pushButtonLoad.clicked.connect(self.loadFile)


        self._ui.scintillaScript = Qsci.QsciScintilla()
        font = QtGui.QFont()
        font.setFamily('Courier')
        font.setPointSize(10)
        self._ui.scintillaScript.setFont(font)
        self._ui.scintillaScript.setMarginsFont(font)
        fontmetrics = QtGui.QFontMetrics(font)
        self._ui.scintillaScript.setMarginWidth(0, fontmetrics.width("000"))
        self._ui.scintillaScript.setMarginLineNumbers(0, True)
        self._ui.scintillaScript.setBraceMatching(Qsci.QsciScintilla.SloppyBraceMatch)
        self._ui.scintillaScript.setCaretLineVisible(True)
        lexer = Qsci.QsciLexerPython()
        lexer.setDefaultFont(font)
        self._ui.scintillaScript.setLexer(lexer)
        self._ui.scintillaScript.SendScintilla(Qsci.QsciScintilla.SCI_STYLESETFONT,1,
                'Courier')
        self._ui.scintillaScript.SendScintilla(Qsci.QsciScintilla.SCI_SETHSCROLLBAR,0)
        self._ui.scintillaScript.setTabWidth(4)
        self._ui.scintillaScript.setIndentationsUseTabs(False)
        self._ui.scintillaScript.setEolMode(Qsci.QsciScintilla.EolMode(Qsci.QsciScintilla.EolUnix))
        self._ui.gridLayout.addWidget(self._ui.scintillaScript, 1, 1)

        script_path = os.path.join(os.path.dirname(__file__), 'scripts')
        for element in os.listdir(script_path):
            element_path = os.path.join(script_path, element)
            if os.path.isfile(os.path.join(element_path)):
                if os.path.splitext(element_path)[1] == '.py':
                    self._add_filename(element_path)
                    

    def _drag_enter_event(self, event):
        if event.mimeData().hasUrls():
            event.accept()
        else:
            event.ignore()

    def _drag_move_event(self, event):
        if event.mimeData().hasUrls():
            event.accept()
        else:
            event.ignore()

    def _drop_event(self, event):
        for url in event.mimeData().urls():
            path = str(url.toLocalFile())
            self._add_filename(path)


    def new(self):
        selected = self._ui.listWidgetFiles.selectedItems()
        if len(selected) > 0:
            self._ui.listWidgetFiles.setItemSelected(selected[0], False)
        self._ui.scintillaScript.clear()


    def loadFile(self):
        filename = str(QtGui.QFileDialog.getOpenFileName(self._window,
            'Open file', '', '*.py'))
        if filename:
            self._add_filename(filename)


    def _add_filename(self, filename):
        existing = self._ui.listWidgetFiles.findItems(filename,
                QtCore.Qt.MatchExactly)
        if len(existing) == 0:
            item = QtGui.QListWidgetItem(filename)
            self._ui.listWidgetFiles.addItem(item)
        else:
            item = existing[0]
        self._ui.listWidgetFiles.setItemSelected(item, True)
        self._update_script(filename)


    def _on_save_clicked(self):
        script = str(self._ui.scintillaScript.text())
        selected = self._ui.listWidgetFiles.selectedItems()
        if len(selected) > 0:
            self._save(str(selected[0].text()), script)
        else:
            filename = str(QtGui.QFileDialog.getSaveFileName(self._window,
                                                             'Save to file',
                                                             '', '*.py'))
            if filename:
                self._save(filename, script)
                self._add_filename(filename)
                

    def _save(self, filename, script):
        f = open(filename, 'w')
        f.write(script)
        f.close()


    def _on_file_clicked(self, item):
        path = str(item.text())
        self._update_script(path)


    def _update_script(self, path):
        f = open(path)
        script = f.read()
        f.close()
        self._ui.scintillaScript.setText(script)


    def run(self):
        if self._ui.pushButtonRun.text() == 'Run':
            self._ui.pushButtonRun.setText('Stop')
            print 'script: Running script.'
            #self._ui.pushButtonRun.setEnabled(False)
            script = str(self._ui.scintillaScript.text())
            self._runner = Runner(self._control, script)
            self._window.connect(self._runner, QtCore.SIGNAL('scriptDone'), self._on_script_done)
            self._runner.start()
        else:
            self._ui.pushButtonRun.setEnabled(False)
            self._ui.pushButtonRun.setText('Stopping')
            self._control.stop = True

    def _on_script_done(self):
        self._runner = None
        self._ui.pushButtonRun.setText('Run')
        self._ui.pushButtonRun.setEnabled(True)
        self._control.stop = False


    def shutDown(self):
        if self._runner:
            print 'Waiting for script to finish.'
            self._runner.wait()



class Runner(QtCore.QThread):
    def __init__(self, control, script):
        QtCore.QThread.__init__(self)
        self.control = control
        self.script = script


    def run(self):
        self.control.run(self.script)
        #self.control.run_interactive(self.script)
        self.emit(QtCore.SIGNAL('scriptDone'))



class MyTextEdit(QtGui.QTextEdit):

    def __init__(self, parent=None):
        QtGui.QTextEdit.__init__(self, parent)

        style = 'font: 9pt "Lucida Console";'
        self.setStyleSheet(style)


    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Tab:
            self.insertPlainText('    ')
        else:
            QtGui.QTextEdit.keyPressEvent(self, event)
