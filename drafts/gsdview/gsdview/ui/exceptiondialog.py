# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'gsdview/ui/exceptiondialog.ui'
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

class Ui_crashDialog(object):
    def setupUi(self, crashDialog):
        crashDialog.setObjectName(_fromUtf8("crashDialog"))
        crashDialog.resize(557, 423)
        crashDialog.setWindowTitle(QtGui.QApplication.translate("crashDialog", "Critical Error", None, QtGui.QApplication.UnicodeUTF8))
        crashDialog.setModal(True)
        self.mainVerticalLayout = QtGui.QVBoxLayout(crashDialog)
        self.mainVerticalLayout.setObjectName(_fromUtf8("mainVerticalLayout"))
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.iconVerticalLayout = QtGui.QVBoxLayout()
        self.iconVerticalLayout.setObjectName(_fromUtf8("iconVerticalLayout"))
        self.iconLabel = QtGui.QLabel(crashDialog)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Maximum, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.iconLabel.sizePolicy().hasHeightForWidth())
        self.iconLabel.setSizePolicy(sizePolicy)
        self.iconLabel.setText(QtGui.QApplication.translate("crashDialog", "Icon", None, QtGui.QApplication.UnicodeUTF8))
        self.iconLabel.setMargin(8)
        self.iconLabel.setObjectName(_fromUtf8("iconLabel"))
        self.iconVerticalLayout.addWidget(self.iconLabel)
        spacerItem = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.iconVerticalLayout.addItem(spacerItem)
        self.horizontalLayout.addLayout(self.iconVerticalLayout)
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.errorMsgLabel = QtGui.QLabel(crashDialog)
        self.errorMsgLabel.setMinimumSize(QtCore.QSize(0, 0))
        self.errorMsgLabel.setText(QtGui.QApplication.translate("crashDialog", "Error message", None, QtGui.QApplication.UnicodeUTF8))
        self.errorMsgLabel.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignTop)
        self.errorMsgLabel.setWordWrap(True)
        self.errorMsgLabel.setOpenExternalLinks(True)
        self.errorMsgLabel.setTextInteractionFlags(QtCore.Qt.LinksAccessibleByMouse|QtCore.Qt.TextSelectableByMouse)
        self.errorMsgLabel.setObjectName(_fromUtf8("errorMsgLabel"))
        self.verticalLayout.addWidget(self.errorMsgLabel)
        self.hline1 = QtGui.QFrame(crashDialog)
        self.hline1.setFrameShape(QtGui.QFrame.HLine)
        self.hline1.setFrameShadow(QtGui.QFrame.Sunken)
        self.hline1.setObjectName(_fromUtf8("hline1"))
        self.verticalLayout.addWidget(self.hline1)
        self.textLabel = QtGui.QLabel(crashDialog)
        self.textLabel.setText(QtGui.QApplication.translate("crashDialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Sans\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Please file a bug report at <a href=\"http://projecturl.example.org\"><span style=\" text-decoration: underline; color:#0000ff;\">projext url</span></a> or report the problem via email to <a href=\"mailto:author_email?subject=[project] Bug report\"><span style=\" text-decoration: underline; color:#0000ff;\">Author Name</span></a>.</p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.textLabel.setWordWrap(True)
        self.textLabel.setTextInteractionFlags(QtCore.Qt.LinksAccessibleByMouse|QtCore.Qt.TextSelectableByMouse)
        self.textLabel.setObjectName(_fromUtf8("textLabel"))
        self.verticalLayout.addWidget(self.textLabel)
        self.hline2 = QtGui.QFrame(crashDialog)
        self.hline2.setFrameShape(QtGui.QFrame.HLine)
        self.hline2.setFrameShadow(QtGui.QFrame.Sunken)
        self.hline2.setObjectName(_fromUtf8("hline2"))
        self.verticalLayout.addWidget(self.hline2)
        self.timestampLabel = QtGui.QLabel(crashDialog)
        self.timestampLabel.setText(QtGui.QApplication.translate("crashDialog", "Timsetamp:", None, QtGui.QApplication.UnicodeUTF8))
        self.timestampLabel.setTextInteractionFlags(QtCore.Qt.LinksAccessibleByMouse|QtCore.Qt.TextSelectableByMouse)
        self.timestampLabel.setObjectName(_fromUtf8("timestampLabel"))
        self.verticalLayout.addWidget(self.timestampLabel)
        self.tracebackGroupBox = QtGui.QGroupBox(crashDialog)
        self.tracebackGroupBox.setTitle(QtGui.QApplication.translate("crashDialog", "Traceback", None, QtGui.QApplication.UnicodeUTF8))
        self.tracebackGroupBox.setCheckable(True)
        self.tracebackGroupBox.setObjectName(_fromUtf8("tracebackGroupBox"))
        self.groupboxVerticalLayout = QtGui.QVBoxLayout(self.tracebackGroupBox)
        self.groupboxVerticalLayout.setObjectName(_fromUtf8("groupboxVerticalLayout"))
        self.tracebackTextEdit = QtGui.QPlainTextEdit(self.tracebackGroupBox)
        self.tracebackTextEdit.setLineWrapMode(QtGui.QPlainTextEdit.NoWrap)
        self.tracebackTextEdit.setReadOnly(True)
        self.tracebackTextEdit.setPlainText(_fromUtf8(""))
        self.tracebackTextEdit.setTextInteractionFlags(QtCore.Qt.TextSelectableByKeyboard|QtCore.Qt.TextSelectableByMouse)
        self.tracebackTextEdit.setObjectName(_fromUtf8("tracebackTextEdit"))
        self.groupboxVerticalLayout.addWidget(self.tracebackTextEdit)
        self.verticalLayout.addWidget(self.tracebackGroupBox)
        spacerItem1 = QtGui.QSpacerItem(20, 0, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem1)
        self.verticalLayout.setStretch(5, 1)
        self.horizontalLayout.addLayout(self.verticalLayout)
        self.horizontalLayout.setStretch(1, 1)
        self.mainVerticalLayout.addLayout(self.horizontalLayout)
        self.buttonBox = QtGui.QDialogButtonBox(crashDialog)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Close|QtGui.QDialogButtonBox.Ignore)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.mainVerticalLayout.addWidget(self.buttonBox)

        self.retranslateUi(crashDialog)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("accepted()")), crashDialog.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("rejected()")), crashDialog.reject)
        QtCore.QObject.connect(self.tracebackGroupBox, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.tracebackTextEdit.setVisible)
        QtCore.QMetaObject.connectSlotsByName(crashDialog)

    def retranslateUi(self, crashDialog):
        pass


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    crashDialog = QtGui.QDialog()
    ui = Ui_crashDialog()
    ui.setupUi(crashDialog)
    crashDialog.show()
    sys.exit(app.exec_())

