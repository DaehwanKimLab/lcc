import time
from datetime import date, datetime

def printTimestamp():
    now = datetime.now()
    nowFormatted = now.strftime("%H:%M:%S")
    print(f"Time: {nowFormatted}")
    return

def printBlockMessage(message, ts:bool = False):
    """ts == timestamp"""
    print("############################################")
    print(message)
    if ts:
        printTimestamp()
    print("############################################")
    return