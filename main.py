import sys
from fbs_runtime.application_context import ApplicationContext
from PyQt5.QtWidgets import QDialog

from DarcyWeisbachFrictionFactor import DarcyWeisbachFrictionFactor, TransitionFlowSolveMethod, TurbulentFlowSolveMethod
from DwFriction import Ui_Dw_Dialog


class DwFrictionApp(QDialog, Ui_Dw_Dialog):
    def __init__(self):
        super(DwFrictionApp, self).__init__()
        self.setupUi(self)
        self.btn_ok.clicked.connect(self.calculate)
        self.setFixedSize(290, 295)

    def calculate(self):
        q = self.m_q.value()
        d = self.m_d.value()
        v = self.m_v.value() * 1e-6
        e = self.m_e.value() * 1e-3
        transition_index = TransitionFlowSolveMethod(self.cb_transition.currentIndex()+1)
        turbulent_index = TurbulentFlowSolveMethod(self.cb_turbulent.currentIndex()+1)
        dw = DarcyWeisbachFrictionFactor(q, d, e=e, v=v, solve_transition_flow_for=transition_index,
                                         solve_turbulent_flow_for=turbulent_index)
        self.te_result.append('<b>Reynolds sayısı, Re = {:.4f}</b>'.format(dw.Re))
        self.te_result.append('<b>Sürtünme katsayısı, f = {:.10f}</b>'.format(dw.f))


class AppContext(ApplicationContext):           # 1. Subclass ApplicationContext
    def run(self):                              # 2. Implement run()
        window = DwFrictionApp()
        # version = self.build_settings['version']
        window.show()
        return self.app.exec_()                 # 3. End run() with this line


if __name__ == '__main__':
    appctxt = AppContext()                      # 4. Instantiate the subclass
    exit_code = appctxt.run()                   # 5. Invoke run()
    sys.exit(exit_code)




