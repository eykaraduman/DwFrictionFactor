import sys

from DarcyWeisbachFrictionFactor import DarcyWeisbachFrictionFactor, TransitionFlowSolveMethod, \
    TurbulentFlowSolveMethod
from DwFriction import Ui_Dw_Dialog
from PyQt5 import QtCore
from PyQt5 import QtWidgets


class DwFrictionApp(QtWidgets.QDialog, Ui_Dw_Dialog):
    def __init__(self):
        super(DwFrictionApp, self).__init__()
        self.setupUi(self)
        # QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create('Cleanlooks'))
        self.btn_ok.clicked.connect(self.calculate)
        self.setFixedSize(290, 295)

    def calculate(self):
        q = self.m_q.value()
        d = self.m_d.value()
        v = self.m_v.value() * 1e-6
        e = self.m_e.value() * 1e-3
        transition_index = TransitionFlowSolveMethod(self.cb_transition.currentIndex()+1)
        turbulant_index = TurbulentFlowSolveMethod(self.cb_turbulent.currentIndex()+1)
        dw = DarcyWeisbachFrictionFactor(q, d, e=e, v=v, solve_transition_flow_for=transition_index,
                                         solve_turbulent_flow_for=turbulant_index)
        self.te_result.append('<b>Reynolds sayısı, Re = {:.4f}</b>'.format(dw.Re))
        self.te_result.append('<b>Sürtünme katsayısı, f = {:.10f}</b>'.format(dw.f))

def main():
    app = QtWidgets.QApplication(sys.argv)
    form = DwFrictionApp()
    form.show()
    app.exec_()

if __name__ == '__main__':
    main()