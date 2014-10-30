# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'gsdview/ui/preferences.ui'
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

class Ui_preferencesDialog(object):
    def setupUi(self, preferencesDialog):
        preferencesDialog.setObjectName(_fromUtf8("preferencesDialog"))
        preferencesDialog.resize(712, 434)
        preferencesDialog.setWindowTitle(QtGui.QApplication.translate("preferencesDialog", "Preferences", None, QtGui.QApplication.UnicodeUTF8))
        self.verticalLayout = QtGui.QVBoxLayout(preferencesDialog)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.listWidget = QtGui.QListWidget(preferencesDialog)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.listWidget.sizePolicy().hasHeightForWidth())
        self.listWidget.setSizePolicy(sizePolicy)
        self.listWidget.setMaximumSize(QtCore.QSize(112, 16777215))
        self.listWidget.setToolTip(_fromUtf8(""))
        self.listWidget.setIconSize(QtCore.QSize(64, 64))
        self.listWidget.setVerticalScrollMode(QtGui.QAbstractItemView.ScrollPerPixel)
        self.listWidget.setMovement(QtGui.QListView.Static)
        self.listWidget.setFlow(QtGui.QListView.TopToBottom)
        self.listWidget.setProperty("isWrapping", False)
        self.listWidget.setResizeMode(QtGui.QListView.Adjust)
        self.listWidget.setSpacing(12)
        self.listWidget.setViewMode(QtGui.QListView.IconMode)
        self.listWidget.setUniformItemSizes(True)
        self.listWidget.setObjectName(_fromUtf8("listWidget"))
        self.horizontalLayout.addWidget(self.listWidget)
        self.stackedWidget = QtGui.QStackedWidget(preferencesDialog)
        self.stackedWidget.setObjectName(_fromUtf8("stackedWidget"))
        self.generalPage = QtGui.QWidget()
        self.generalPage.setObjectName(_fromUtf8("generalPage"))
        self.generalVerticalLayout = QtGui.QVBoxLayout(self.generalPage)
        self.generalVerticalLayout.setSpacing(3)
        self.generalVerticalLayout.setObjectName(_fromUtf8("generalVerticalLayout"))
        self.stackedWidget.addWidget(self.generalPage)
        self.horizontalLayout.addWidget(self.stackedWidget)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.buttonBox = QtGui.QDialogButtonBox(preferencesDialog)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Apply|QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.verticalLayout.addWidget(self.buttonBox)

        self.retranslateUi(preferencesDialog)
        self.stackedWidget.setCurrentIndex(0)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("accepted()")), preferencesDialog.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("rejected()")), preferencesDialog.reject)
        QtCore.QMetaObject.connectSlotsByName(preferencesDialog)

    def retranslateUi(self, preferencesDialog):
        pass


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    preferencesDialog = QtGui.QDialog()
    ui = Ui_preferencesDialog()
    ui.setupUi(preferencesDialog)
    preferencesDialog.show()
    sys.exit(app.exec_())

