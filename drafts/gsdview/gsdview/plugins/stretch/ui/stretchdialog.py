# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'gsdview/plugins/stretch/ui/stretchdialog.ui'
#
# Created: Sun Nov 13 21:17:58 2011
#      by: PyQt4 UI code generator 4.8.5
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_stretchingDialog(object):
    def setupUi(self, stretchingDialog):
        stretchingDialog.setObjectName(_fromUtf8("stretchingDialog"))
        stretchingDialog.resize(617, 100)
        stretchingDialog.setWindowTitle(QtGui.QApplication.translate("stretchingDialog", "Linear stretching", None, QtGui.QApplication.UnicodeUTF8))
        self.mainLayout = QtGui.QVBoxLayout(stretchingDialog)
        self.mainLayout.setObjectName(_fromUtf8("mainLayout"))
        spacerItem = QtGui.QSpacerItem(20, 53, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.mainLayout.addItem(spacerItem)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.checkBox = QtGui.QCheckBox(stretchingDialog)
        self.checkBox.setText(QtGui.QApplication.translate("stretchingDialog", "Advanced", None, QtGui.QApplication.UnicodeUTF8))
        self.checkBox.setObjectName(_fromUtf8("checkBox"))
        self.horizontalLayout.addWidget(self.checkBox)
        spacerItem1 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem1)
        self.resetButton = QtGui.QPushButton(stretchingDialog)
        self.resetButton.setText(QtGui.QApplication.translate("stretchingDialog", "Reset", None, QtGui.QApplication.UnicodeUTF8))
        self.resetButton.setObjectName(_fromUtf8("resetButton"))
        self.horizontalLayout.addWidget(self.resetButton)
        self.mainLayout.addLayout(self.horizontalLayout)

        self.retranslateUi(stretchingDialog)
        QtCore.QMetaObject.connectSlotsByName(stretchingDialog)

    def retranslateUi(self, stretchingDialog):
        pass


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    stretchingDialog = QtGui.QDialog()
    ui = Ui_stretchingDialog()
    ui.setupUi(stretchingDialog)
    stretchingDialog.show()
    sys.exit(app.exec_())

