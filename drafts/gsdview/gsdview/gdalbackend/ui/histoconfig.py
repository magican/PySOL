# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'gsdview/gdalbackend/ui/histoconfig.ui'
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

class Ui_histogramConfigDialog(object):
    def setupUi(self, histogramConfigDialog):
        histogramConfigDialog.setObjectName(_fromUtf8("histogramConfigDialog"))
        histogramConfigDialog.resize(350, 208)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(histogramConfigDialog.sizePolicy().hasHeightForWidth())
        histogramConfigDialog.setSizePolicy(sizePolicy)
        histogramConfigDialog.setWindowTitle(QtGui.QApplication.translate("histogramConfigDialog", "Histogram Configuration", None, QtGui.QApplication.UnicodeUTF8))
        histogramConfigDialog.setModal(True)
        self.verticalLayout = QtGui.QVBoxLayout(histogramConfigDialog)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.formLayout = QtGui.QFormLayout()
        self.formLayout.setObjectName(_fromUtf8("formLayout"))
        self.minimumLabel = QtGui.QLabel(histogramConfigDialog)
        self.minimumLabel.setText(QtGui.QApplication.translate("histogramConfigDialog", "Minimum:", None, QtGui.QApplication.UnicodeUTF8))
        self.minimumLabel.setObjectName(_fromUtf8("minimumLabel"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.LabelRole, self.minimumLabel)
        self.minSpinBox = QtGui.QDoubleSpinBox(histogramConfigDialog)
        self.minSpinBox.setMinimum(-32768.0)
        self.minSpinBox.setMaximum(65536.0)
        self.minSpinBox.setObjectName(_fromUtf8("minSpinBox"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.FieldRole, self.minSpinBox)
        self.maximumLabel = QtGui.QLabel(histogramConfigDialog)
        self.maximumLabel.setText(QtGui.QApplication.translate("histogramConfigDialog", "Maximum:", None, QtGui.QApplication.UnicodeUTF8))
        self.maximumLabel.setObjectName(_fromUtf8("maximumLabel"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.LabelRole, self.maximumLabel)
        self.maxSpinBox = QtGui.QDoubleSpinBox(histogramConfigDialog)
        self.maxSpinBox.setMaximum(65536.0)
        self.maxSpinBox.setProperty("value", 256.0)
        self.maxSpinBox.setObjectName(_fromUtf8("maxSpinBox"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.FieldRole, self.maxSpinBox)
        self.nBucketsLabel = QtGui.QLabel(histogramConfigDialog)
        self.nBucketsLabel.setText(QtGui.QApplication.translate("histogramConfigDialog", "Nubber of buckets:", None, QtGui.QApplication.UnicodeUTF8))
        self.nBucketsLabel.setObjectName(_fromUtf8("nBucketsLabel"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.LabelRole, self.nBucketsLabel)
        self.nBucketsSpinBox = QtGui.QSpinBox(histogramConfigDialog)
        self.nBucketsSpinBox.setToolTip(QtGui.QApplication.translate("histogramConfigDialog", "The number of buckets in the histogram.", None, QtGui.QApplication.UnicodeUTF8))
        self.nBucketsSpinBox.setMinimum(4)
        self.nBucketsSpinBox.setMaximum(10000)
        self.nBucketsSpinBox.setProperty("value", 256)
        self.nBucketsSpinBox.setObjectName(_fromUtf8("nBucketsSpinBox"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.FieldRole, self.nBucketsSpinBox)
        self.outOfRangeCheckBox = QtGui.QCheckBox(histogramConfigDialog)
        self.outOfRangeCheckBox.setToolTip(QtGui.QApplication.translate("histogramConfigDialog", "If checked values below the histogram range will mapped into the first value,\n"
"and values above will be mapped into the last one otherwise out of range \n"
"values are discarded.", None, QtGui.QApplication.UnicodeUTF8))
        self.outOfRangeCheckBox.setText(QtGui.QApplication.translate("histogramConfigDialog", "Include out of range", None, QtGui.QApplication.UnicodeUTF8))
        self.outOfRangeCheckBox.setObjectName(_fromUtf8("outOfRangeCheckBox"))
        self.formLayout.setWidget(3, QtGui.QFormLayout.FieldRole, self.outOfRangeCheckBox)
        self.approxCheckBox = QtGui.QCheckBox(histogramConfigDialog)
        self.approxCheckBox.setEnabled(False)
        self.approxCheckBox.setToolTip(QtGui.QApplication.translate("histogramConfigDialog", "If checked allows approximate or incomplete \n"
"histogram computation (faster).", None, QtGui.QApplication.UnicodeUTF8))
        self.approxCheckBox.setText(QtGui.QApplication.translate("histogramConfigDialog", "Approx OK", None, QtGui.QApplication.UnicodeUTF8))
        self.approxCheckBox.setChecked(True)
        self.approxCheckBox.setObjectName(_fromUtf8("approxCheckBox"))
        self.formLayout.setWidget(4, QtGui.QFormLayout.FieldRole, self.approxCheckBox)
        self.verticalLayout.addLayout(self.formLayout)
        self.dialogButtonBox = QtGui.QDialogButtonBox(histogramConfigDialog)
        self.dialogButtonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.dialogButtonBox.setObjectName(_fromUtf8("dialogButtonBox"))
        self.verticalLayout.addWidget(self.dialogButtonBox)

        self.retranslateUi(histogramConfigDialog)
        QtCore.QObject.connect(self.dialogButtonBox, QtCore.SIGNAL(_fromUtf8("accepted()")), histogramConfigDialog.accept)
        QtCore.QObject.connect(self.dialogButtonBox, QtCore.SIGNAL(_fromUtf8("rejected()")), histogramConfigDialog.reject)
        QtCore.QMetaObject.connectSlotsByName(histogramConfigDialog)

    def retranslateUi(self, histogramConfigDialog):
        pass


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    histogramConfigDialog = QtGui.QDialog()
    ui = Ui_histogramConfigDialog()
    ui.setupUi(histogramConfigDialog)
    histogramConfigDialog.show()
    sys.exit(app.exec_())

