# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'gsdview/plugins/stretch/ui/doubleslider.ui'
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

class Ui_doubleslider(object):
    def setupUi(self, doubleslider):
        doubleslider.setObjectName(_fromUtf8("doubleslider"))
        doubleslider.resize(617, 134)
        doubleslider.setWindowTitle(QtGui.QApplication.translate("doubleslider", "Form", None, QtGui.QApplication.UnicodeUTF8))
        self.verticalLayout = QtGui.QVBoxLayout(doubleslider)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.lowLayout = QtGui.QHBoxLayout()
        self.lowLayout.setObjectName(_fromUtf8("lowLayout"))
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.lowLayout.addItem(spacerItem)
        self.lowSpinBox = QtGui.QDoubleSpinBox(doubleslider)
        self.lowSpinBox.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.lowSpinBox.setMinimum(-100000.0)
        self.lowSpinBox.setMaximum(100000.0)
        self.lowSpinBox.setObjectName(_fromUtf8("lowSpinBox"))
        self.lowLayout.addWidget(self.lowSpinBox)
        spacerItem1 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.lowLayout.addItem(spacerItem1)
        self.verticalLayout.addLayout(self.lowLayout)
        self.stretchLayout = QtGui.QHBoxLayout()
        self.stretchLayout.setObjectName(_fromUtf8("stretchLayout"))
        self.minSpinBox = QtGui.QDoubleSpinBox(doubleslider)
        self.minSpinBox.setMaximumSize(QtCore.QSize(100, 16777215))
        self.minSpinBox.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.minSpinBox.setMinimum(-100000.0)
        self.minSpinBox.setMaximum(100000.0)
        self.minSpinBox.setProperty("value", 0.0)
        self.minSpinBox.setObjectName(_fromUtf8("minSpinBox"))
        self.stretchLayout.addWidget(self.minSpinBox)
        self.slidersLayout = QtGui.QVBoxLayout()
        self.slidersLayout.setObjectName(_fromUtf8("slidersLayout"))
        self.lowSlider = QtGui.QSlider(doubleslider)
        self.lowSlider.setMaximum(255)
        self.lowSlider.setProperty("value", 0)
        self.lowSlider.setOrientation(QtCore.Qt.Horizontal)
        self.lowSlider.setObjectName(_fromUtf8("lowSlider"))
        self.slidersLayout.addWidget(self.lowSlider)
        self.highSlider = QtGui.QSlider(doubleslider)
        self.highSlider.setMaximum(255)
        self.highSlider.setProperty("value", 100)
        self.highSlider.setOrientation(QtCore.Qt.Horizontal)
        self.highSlider.setObjectName(_fromUtf8("highSlider"))
        self.slidersLayout.addWidget(self.highSlider)
        self.stretchLayout.addLayout(self.slidersLayout)
        self.maxSpinBox = QtGui.QDoubleSpinBox(doubleslider)
        self.maxSpinBox.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.maxSpinBox.setMaximum(100000.0)
        self.maxSpinBox.setProperty("value", 255.0)
        self.maxSpinBox.setObjectName(_fromUtf8("maxSpinBox"))
        self.stretchLayout.addWidget(self.maxSpinBox)
        self.verticalLayout.addLayout(self.stretchLayout)
        self.highLayout = QtGui.QHBoxLayout()
        self.highLayout.setObjectName(_fromUtf8("highLayout"))
        spacerItem2 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.highLayout.addItem(spacerItem2)
        self.highSpinBox = QtGui.QDoubleSpinBox(doubleslider)
        self.highSpinBox.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.highSpinBox.setMaximum(100000.0)
        self.highSpinBox.setProperty("value", 100.0)
        self.highSpinBox.setObjectName(_fromUtf8("highSpinBox"))
        self.highLayout.addWidget(self.highSpinBox)
        spacerItem3 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.highLayout.addItem(spacerItem3)
        self.verticalLayout.addLayout(self.highLayout)

        self.retranslateUi(doubleslider)
        QtCore.QMetaObject.connectSlotsByName(doubleslider)

    def retranslateUi(self, doubleslider):
        pass


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    doubleslider = QtGui.QWidget()
    ui = Ui_doubleslider()
    ui.setupUi(doubleslider)
    doubleslider.show()
    sys.exit(app.exec_())

