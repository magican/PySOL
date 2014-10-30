# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'gsdview/ui/general-page.ui'
#
# Created: Sun Nov 13 21:17:57 2011
#      by: PyQt4 UI code generator 4.8.5
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_generalPreferencesWidget(object):
    def setupUi(self, generalPreferencesWidget):
        generalPreferencesWidget.setObjectName(_fromUtf8("generalPreferencesWidget"))
        generalPreferencesWidget.resize(591, 326)
        generalPreferencesWidget.setWindowTitle(QtGui.QApplication.translate("generalPreferencesWidget", "Form", None, QtGui.QApplication.UnicodeUTF8))
        self.verticalLayout = QtGui.QVBoxLayout(generalPreferencesWidget)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.preferencesGroupBox = QtGui.QGroupBox(generalPreferencesWidget)
        self.preferencesGroupBox.setTitle(QtGui.QApplication.translate("generalPreferencesWidget", "Preferences", None, QtGui.QApplication.UnicodeUTF8))
        self.preferencesGroupBox.setObjectName(_fromUtf8("preferencesGroupBox"))
        self.preferencesGridLayout = QtGui.QGridLayout(self.preferencesGroupBox)
        self.preferencesGridLayout.setHorizontalSpacing(12)
        self.preferencesGridLayout.setObjectName(_fromUtf8("preferencesGridLayout"))
        self.loglevelLabel = QtGui.QLabel(self.preferencesGroupBox)
        self.loglevelLabel.setText(QtGui.QApplication.translate("generalPreferencesWidget", "Log level", None, QtGui.QApplication.UnicodeUTF8))
        self.loglevelLabel.setObjectName(_fromUtf8("loglevelLabel"))
        self.preferencesGridLayout.addWidget(self.loglevelLabel, 0, 0, 1, 1)
        self.cachedirLabel = QtGui.QLabel(self.preferencesGroupBox)
        self.cachedirLabel.setText(QtGui.QApplication.translate("generalPreferencesWidget", "Cache directory location", None, QtGui.QApplication.UnicodeUTF8))
        self.cachedirLabel.setObjectName(_fromUtf8("cachedirLabel"))
        self.preferencesGridLayout.addWidget(self.cachedirLabel, 1, 0, 1, 1)
        self.loglevelHorizontalLayout = QtGui.QHBoxLayout()
        self.loglevelHorizontalLayout.setObjectName(_fromUtf8("loglevelHorizontalLayout"))
        self.loglevelComboBox = QtGui.QComboBox(self.preferencesGroupBox)
        self.loglevelComboBox.setObjectName(_fromUtf8("loglevelComboBox"))
        self.loglevelComboBox.addItem(_fromUtf8(""))
        self.loglevelComboBox.setItemText(0, QtGui.QApplication.translate("generalPreferencesWidget", "DEBUG", None, QtGui.QApplication.UnicodeUTF8))
        self.loglevelComboBox.addItem(_fromUtf8(""))
        self.loglevelComboBox.setItemText(1, QtGui.QApplication.translate("generalPreferencesWidget", "INFO", None, QtGui.QApplication.UnicodeUTF8))
        self.loglevelComboBox.addItem(_fromUtf8(""))
        self.loglevelComboBox.setItemText(2, QtGui.QApplication.translate("generalPreferencesWidget", "WARNING", None, QtGui.QApplication.UnicodeUTF8))
        self.loglevelComboBox.addItem(_fromUtf8(""))
        self.loglevelComboBox.setItemText(3, QtGui.QApplication.translate("generalPreferencesWidget", "ERROR", None, QtGui.QApplication.UnicodeUTF8))
        self.loglevelComboBox.addItem(_fromUtf8(""))
        self.loglevelComboBox.setItemText(4, QtGui.QApplication.translate("generalPreferencesWidget", "CRITICAL", None, QtGui.QApplication.UnicodeUTF8))
        self.loglevelHorizontalLayout.addWidget(self.loglevelComboBox)
        spacerItem = QtGui.QSpacerItem(138, 24, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.loglevelHorizontalLayout.addItem(spacerItem)
        self.preferencesGridLayout.addLayout(self.loglevelHorizontalLayout, 0, 1, 1, 1)
        self.verticalLayout.addWidget(self.preferencesGroupBox)
        self.filedialogGroupBox = QtGui.QGroupBox(generalPreferencesWidget)
        self.filedialogGroupBox.setTitle(QtGui.QApplication.translate("generalPreferencesWidget", "File DIalog", None, QtGui.QApplication.UnicodeUTF8))
        self.filedialogGroupBox.setObjectName(_fromUtf8("filedialogGroupBox"))
        self.fileDialogHorizontalLayout = QtGui.QHBoxLayout(self.filedialogGroupBox)
        self.fileDialogHorizontalLayout.setObjectName(_fromUtf8("fileDialogHorizontalLayout"))
        self.workdirLabel = QtGui.QLabel(self.filedialogGroupBox)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.workdirLabel.sizePolicy().hasHeightForWidth())
        self.workdirLabel.setSizePolicy(sizePolicy)
        self.workdirLabel.setText(QtGui.QApplication.translate("generalPreferencesWidget", "Working directory", None, QtGui.QApplication.UnicodeUTF8))
        self.workdirLabel.setObjectName(_fromUtf8("workdirLabel"))
        self.fileDialogHorizontalLayout.addWidget(self.workdirLabel)
        self.verticalLayout.addWidget(self.filedialogGroupBox)
        spacerItem1 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem1)

        self.retranslateUi(generalPreferencesWidget)
        QtCore.QMetaObject.connectSlotsByName(generalPreferencesWidget)

    def retranslateUi(self, generalPreferencesWidget):
        pass


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    generalPreferencesWidget = QtGui.QWidget()
    ui = Ui_generalPreferencesWidget()
    ui.setupUi(generalPreferencesWidget)
    generalPreferencesWidget.show()
    sys.exit(app.exec_())

