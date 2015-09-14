__author__ = 'juncheng'




class SPEMD:
    def __init__(self):
        self.name="SPEMD"
        self.critRad=0
        self.ellipticity = 0
        self.axisRatio = 0
        self.pa = 0
        self.gamma = 0
        self.coreRad = 0
        self.centerX=0
        self.centerY=0

class SIP:
    def __init__(self, critRad, axisRatio, pa, centerX, centerY ):
        self.name="SIE"
        self.critRad=critRad
        self.axisRatio = axisRatio
        self.pa = pa
        self.centerX = centerX
        self.centerY = centerY


def main():
    return "Nothing"


if __name__ == "__main__":
    main()