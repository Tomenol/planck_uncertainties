# -*- coding: utf-8 -*-
"""
Created on Tue May 25 10:55:43 2021

@author: Thomas Maynadié
"""

from datetime import datetime
from enum import Enum

DEBUG_MODE = False

class EVSErrCode(Enum):
    EVS_ERROR = 1
    EVS_STATUS = 0

def debug(msg, errorCode=EVSErrCode.EVS_STATUS):
    HEADER = datetime.now().strftime("%H:%M:%S.%f")[:-3] + " >> "
    
    if errorCode == EVSErrCode.EVS_ERROR: 
        ERROR_MSG = "[EVS ERROR] "        
    if errorCode == EVSErrCode.EVS_STATUS: 
        ERROR_MSG = "[EVS STATUS] "
    
    if DEBUG_MODE is True: print(HEADER + ERROR_MSG + msg)