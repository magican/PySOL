# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'gsdview/gdalbackend/ui/metadata.ui'
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

class Ui_metadataForm(object):
    def setupUi(self, metadataForm):
        metadataForm.setObjectName(_fromUtf8("metadataForm"))
        metadataForm.resize(424, 355)
        metadataForm.setWindowTitle(QtGui.QApplication.translate("metadataForm", "Form", None, QtGui.QApplication.UnicodeUTF8))
        self.verticalLayout = QtGui.QVBoxLayout(metadataForm)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.metadataHorizontalLayout = QtGui.QHBoxLayout()
        self.metadataHorizontalLayout.setObjectName(_fromUtf8("metadataHorizontalLayout"))
        self.domainLabel = QtGui.QLabel(metadataForm)
        self.domainLabel.setText(QtGui.QApplication.translate("metadataForm", "Domain:", None, QtGui.QApplication.UnicodeUTF8))
        self.domainLabel.setObjectName(_fromUtf8("domainLabel"))
        self.metadataHorizontalLayout.addWidget(self.domainLabel)
        self.domainComboBox = QtGui.QComboBox(metadataForm)
        self.domainComboBox.setEditable(True)
        self.domainComboBox.setObjectName(_fromUtf8("domainComboBox"))
        self.domainComboBox.addItem(_fromUtf8(""))
        self.domainComboBox.setItemText(0, _fromUtf8(""))
        self.domainComboBox.addItem(_fromUtf8(""))
        self.domainComboBox.setItemText(1, QtGui.QApplication.translate("metadataForm", "SUBDATASETS", None, QtGui.QApplication.UnicodeUTF8))
        self.domainComboBox.addItem(_fromUtf8(""))
        self.domainComboBox.setItemText(2, QtGui.QApplication.translate("metadataForm", "IMAGE_STRUCTURE", None, QtGui.QApplication.UnicodeUTF8))
        self.domainComboBox.addItem(_fromUtf8(""))
        self.domainComboBox.setItemText(3, QtGui.QApplication.translate("metadataForm", "RPC", None, QtGui.QApplication.UnicodeUTF8))
        self.metadataHorizontalLayout.addWidget(self.domainComboBox)
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.metadataHorizontalLayout.addItem(spacerItem)
        self.metadataNumLabel = QtGui.QLabel(metadataForm)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.metadataNumLabel.sizePolicy().hasHeightForWidth())
        self.metadataNumLabel.setSizePolicy(sizePolicy)
        self.metadataNumLabel.setText(QtGui.QApplication.translate("metadataForm", "Number of metadata:", None, QtGui.QApplication.UnicodeUTF8))
        self.metadataNumLabel.setObjectName(_fromUtf8("metadataNumLabel"))
        self.metadataHorizontalLayout.addWidget(self.metadataNumLabel)
        self.metadataNumValue = QtGui.QLabel(metadataForm)
        self.metadataNumValue.setText(QtGui.QApplication.translate("metadataForm", "0", None, QtGui.QApplication.UnicodeUTF8))
        self.metadataNumValue.setTextInteractionFlags(QtCore.Qt.LinksAccessibleByMouse|QtCore.Qt.TextSelectableByKeyboard|QtCore.Qt.TextSelectableByMouse)
        self.metadataNumValue.setObjectName(_fromUtf8("metadataNumValue"))
        self.metadataHorizontalLayout.addWidget(self.metadataNumValue)
        self.verticalLayout.addLayout(self.metadataHorizontalLayout)
        self.tableWidget = QtGui.QTableWidget(metadataForm)
        self.tableWidget.setContextMenuPolicy(QtCore.Qt.ActionsContextMenu)
        self.tableWidget.setEditTriggers(QtGui.QAbstractItemView.NoEditTriggers)
        self.tableWidget.setHorizontalScrollMode(QtGui.QAbstractItemView.ScrollPerPixel)
        self.tableWidget.setColumnCount(2)
        self.tableWidget.setObjectName(_fromUtf8("tableWidget"))
        self.tableWidget.setRowCount(0)
        item = QtGui.QTableWidgetItem()
        item.setText(QtGui.QApplication.translate("metadataForm", "Name", None, QtGui.QApplication.UnicodeUTF8))
        self.tableWidget.setHorizontalHeaderItem(0, item)
        item = QtGui.QTableWidgetItem()
        item.setText(QtGui.QApplication.translate("metadataForm", "Value", None, QtGui.QApplication.UnicodeUTF8))
        self.tableWidget.setHorizontalHeaderItem(1, item)
        self.tableWidget.horizontalHeader().setStretchLastSection(True)
        self.verticalLayout.addWidget(self.tableWidget)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.exportButton = QtGui.QPushButton(metadataForm)
        self.exportButton.setToolTip(QtGui.QApplication.translate("metadataForm", "Save matadata to file.", None, QtGui.QApplication.UnicodeUTF8))
        self.exportButton.setStatusTip(QtGui.QApplication.translate("metadataForm", "Save matadata to file.", None, QtGui.QApplication.UnicodeUTF8))
        self.exportButton.setText(QtGui.QApplication.translate("metadataForm", "Export", None, QtGui.QApplication.UnicodeUTF8))
        self.exportButton.setShortcut(QtGui.QApplication.translate("metadataForm", "Ctrl+S", None, QtGui.QApplication.UnicodeUTF8))
        self.exportButton.setObjectName(_fromUtf8("exportButton"))
        self.horizontalLayout.addWidget(self.exportButton)
        self.printButton = QtGui.QPushButton(metadataForm)
        self.printButton.setStatusTip(QtGui.QApplication.translate("metadataForm", "Print metadata.", None, QtGui.QApplication.UnicodeUTF8))
        self.printButton.setWhatsThis(QtGui.QApplication.translate("metadataForm", "Print metadata.", None, QtGui.QApplication.UnicodeUTF8))
        self.printButton.setText(QtGui.QApplication.translate("metadataForm", "&Print", None, QtGui.QApplication.UnicodeUTF8))
        self.printButton.setShortcut(QtGui.QApplication.translate("metadataForm", "Ctrl+P", None, QtGui.QApplication.UnicodeUTF8))
        self.printButton.setObjectName(_fromUtf8("printButton"))
        self.horizontalLayout.addWidget(self.printButton)
        spacerItem1 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem1)
        self.verticalLayout.addLayout(self.horizontalLayout)

        self.retranslateUi(metadataForm)
        QtCore.QMetaObject.connectSlotsByName(metadataForm)

    def retranslateUi(self, metadataForm):
        self.tableWidget.setSortingEnabled(True)
        item = self.tableWidget.horizontalHeaderItem(0)
        item = self.tableWidget.horizontalHeaderItem(1)


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    metadataForm = QtGui.QWidget()
    ui = Ui_metadataForm()
    ui.setupUi(metadataForm)
    metadataForm.show()
    sys.exit(app.exec_())

