# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'gsdview/gdalbackend/ui/gdalpage.ui'
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

class Ui_gdalPreferencesWidget(object):
    def setupUi(self, gdalPreferencesWidget):
        gdalPreferencesWidget.setObjectName(_fromUtf8("gdalPreferencesWidget"))
        gdalPreferencesWidget.resize(591, 581)
        gdalPreferencesWidget.setWindowTitle(QtGui.QApplication.translate("gdalPreferencesWidget", "Form", None, QtGui.QApplication.UnicodeUTF8))
        self.verticalLayout = QtGui.QVBoxLayout(gdalPreferencesWidget)
        self.verticalLayout.setSpacing(0)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.infoHorizontalLayout = QtGui.QHBoxLayout()
        self.infoHorizontalLayout.setObjectName(_fromUtf8("infoHorizontalLayout"))
        spacerItem = QtGui.QSpacerItem(188, 24, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.infoHorizontalLayout.addItem(spacerItem)
        self.infoButton = QtGui.QPushButton(gdalPreferencesWidget)
        self.infoButton.setText(QtGui.QApplication.translate("gdalPreferencesWidget", "Show detailed GDAL info", None, QtGui.QApplication.UnicodeUTF8))
        self.infoButton.setObjectName(_fromUtf8("infoButton"))
        self.infoHorizontalLayout.addWidget(self.infoButton)
        self.verticalLayout.addLayout(self.infoHorizontalLayout)
        self.optionsGroupBox = QtGui.QGroupBox(gdalPreferencesWidget)
        self.optionsGroupBox.setTitle(QtGui.QApplication.translate("gdalPreferencesWidget", "Configuration options", None, QtGui.QApplication.UnicodeUTF8))
        self.optionsGroupBox.setObjectName(_fromUtf8("optionsGroupBox"))
        self.optionsGridLayout = QtGui.QGridLayout(self.optionsGroupBox)
        self.optionsGridLayout.setVerticalSpacing(2)
        self.optionsGridLayout.setObjectName(_fromUtf8("optionsGridLayout"))
        self.cacheCheckBox = QtGui.QCheckBox(self.optionsGroupBox)
        self.cacheCheckBox.setText(QtGui.QApplication.translate("gdalPreferencesWidget", "GDAL_CACHEMAX", None, QtGui.QApplication.UnicodeUTF8))
        self.cacheCheckBox.setObjectName(_fromUtf8("cacheCheckBox"))
        self.optionsGridLayout.addWidget(self.cacheCheckBox, 0, 0, 1, 1)
        self.cacheHorizontalLayout = QtGui.QHBoxLayout()
        self.cacheHorizontalLayout.setObjectName(_fromUtf8("cacheHorizontalLayout"))
        self.cacheSpinBox = QtGui.QSpinBox(self.optionsGroupBox)
        self.cacheSpinBox.setEnabled(False)
        self.cacheSpinBox.setSuffix(QtGui.QApplication.translate("gdalPreferencesWidget", " MB", None, QtGui.QApplication.UnicodeUTF8))
        self.cacheSpinBox.setMaximum(1000)
        self.cacheSpinBox.setSingleStep(10)
        self.cacheSpinBox.setProperty("value", 40)
        self.cacheSpinBox.setObjectName(_fromUtf8("cacheSpinBox"))
        self.cacheHorizontalLayout.addWidget(self.cacheSpinBox)
        spacerItem1 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.cacheHorizontalLayout.addItem(spacerItem1)
        self.optionsGridLayout.addLayout(self.cacheHorizontalLayout, 0, 1, 1, 1)
        self.gdalDataCheckBox = QtGui.QCheckBox(self.optionsGroupBox)
        self.gdalDataCheckBox.setText(QtGui.QApplication.translate("gdalPreferencesWidget", "GDAL_DATA", None, QtGui.QApplication.UnicodeUTF8))
        self.gdalDataCheckBox.setObjectName(_fromUtf8("gdalDataCheckBox"))
        self.optionsGridLayout.addWidget(self.gdalDataCheckBox, 1, 0, 1, 1)
        self.skipCheckBox = QtGui.QCheckBox(self.optionsGroupBox)
        self.skipCheckBox.setText(QtGui.QApplication.translate("gdalPreferencesWidget", "GDAL_SKIP", None, QtGui.QApplication.UnicodeUTF8))
        self.skipCheckBox.setObjectName(_fromUtf8("skipCheckBox"))
        self.optionsGridLayout.addWidget(self.skipCheckBox, 2, 0, 1, 1)
        self.skipLineEdit = QtGui.QLineEdit(self.optionsGroupBox)
        self.skipLineEdit.setEnabled(False)
        self.skipLineEdit.setObjectName(_fromUtf8("skipLineEdit"))
        self.optionsGridLayout.addWidget(self.skipLineEdit, 2, 1, 1, 1)
        self.gdalDriverPathCheckBox = QtGui.QCheckBox(self.optionsGroupBox)
        self.gdalDriverPathCheckBox.setText(QtGui.QApplication.translate("gdalPreferencesWidget", "GDAL_DRIVER_PATH", None, QtGui.QApplication.UnicodeUTF8))
        self.gdalDriverPathCheckBox.setObjectName(_fromUtf8("gdalDriverPathCheckBox"))
        self.optionsGridLayout.addWidget(self.gdalDriverPathCheckBox, 3, 0, 1, 1)
        self.ogrDriverPathCheckBox = QtGui.QCheckBox(self.optionsGroupBox)
        self.ogrDriverPathCheckBox.setText(QtGui.QApplication.translate("gdalPreferencesWidget", "OGR_DRIVER_PATH", None, QtGui.QApplication.UnicodeUTF8))
        self.ogrDriverPathCheckBox.setObjectName(_fromUtf8("ogrDriverPathCheckBox"))
        self.optionsGridLayout.addWidget(self.ogrDriverPathCheckBox, 4, 0, 1, 1)
        self.verticalLayout.addWidget(self.optionsGroupBox)
        self.extraOptGroupBox = QtGui.QGroupBox(gdalPreferencesWidget)
        self.extraOptGroupBox.setTitle(QtGui.QApplication.translate("gdalPreferencesWidget", "Extra configuration optons", None, QtGui.QApplication.UnicodeUTF8))
        self.extraOptGroupBox.setCheckable(True)
        self.extraOptGroupBox.setChecked(False)
        self.extraOptGroupBox.setObjectName(_fromUtf8("extraOptGroupBox"))
        self.extraOptVerticalLayout = QtGui.QVBoxLayout(self.extraOptGroupBox)
        self.extraOptVerticalLayout.setObjectName(_fromUtf8("extraOptVerticalLayout"))
        self.extraOptTableWidget = QtGui.QTableWidget(self.extraOptGroupBox)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.extraOptTableWidget.sizePolicy().hasHeightForWidth())
        self.extraOptTableWidget.setSizePolicy(sizePolicy)
        self.extraOptTableWidget.setContextMenuPolicy(QtCore.Qt.ActionsContextMenu)
        self.extraOptTableWidget.setVerticalScrollMode(QtGui.QAbstractItemView.ScrollPerPixel)
        self.extraOptTableWidget.setHorizontalScrollMode(QtGui.QAbstractItemView.ScrollPerPixel)
        self.extraOptTableWidget.setObjectName(_fromUtf8("extraOptTableWidget"))
        self.extraOptTableWidget.setColumnCount(2)
        self.extraOptTableWidget.setRowCount(0)
        item = QtGui.QTableWidgetItem()
        item.setText(QtGui.QApplication.translate("gdalPreferencesWidget", "Key", None, QtGui.QApplication.UnicodeUTF8))
        self.extraOptTableWidget.setHorizontalHeaderItem(0, item)
        item = QtGui.QTableWidgetItem()
        item.setText(QtGui.QApplication.translate("gdalPreferencesWidget", "Value", None, QtGui.QApplication.UnicodeUTF8))
        self.extraOptTableWidget.setHorizontalHeaderItem(1, item)
        self.extraOptTableWidget.horizontalHeader().setStretchLastSection(True)
        self.extraOptVerticalLayout.addWidget(self.extraOptTableWidget)
        self.verticalLayout.addWidget(self.extraOptGroupBox)

        self.retranslateUi(gdalPreferencesWidget)
        QtCore.QObject.connect(self.cacheCheckBox, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.cacheSpinBox.setEnabled)
        QtCore.QObject.connect(self.skipCheckBox, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.skipLineEdit.setEnabled)
        QtCore.QMetaObject.connectSlotsByName(gdalPreferencesWidget)

    def retranslateUi(self, gdalPreferencesWidget):
        item = self.extraOptTableWidget.horizontalHeaderItem(0)
        item = self.extraOptTableWidget.horizontalHeaderItem(1)


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    gdalPreferencesWidget = QtGui.QWidget()
    ui = Ui_gdalPreferencesWidget()
    ui.setupUi(gdalPreferencesWidget)
    gdalPreferencesWidget.show()
    sys.exit(app.exec_())

