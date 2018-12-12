# -*- coding: utf-8 -*-

from mdstudio.runner import main
from mdstudio_smartcyp.wamp_services import SmartCypWampApi

# from mdstudio_smartcyp.smartcyp_run import SmartCypRunner
# s = SmartCypRunner()
# d = s.run(mol='O=Cc1ccc(s1)c2cccnc2', is_smiles=True)
# print(d)

if __name__ == '__main__':
    main(SmartCypWampApi)
