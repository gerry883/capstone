# -*- coding: utf-8 -*-
#  Copyright (c) 2016-2017, Zhijiang Yao, Jie Dong and Dongsheng Cao
#  All rights reserved.
#  This file is part of the PyBioMed.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the PyBioMed source tree.
"""
A class used for computing different types of protein descriptors! 

You can freely use and distribute it. If you have any problem, 

you could contact with us timely.

Authors: Zhijiang Yao and Dongsheng Cao.

Date: 2016.06.04

Email: gadsby@163.com

"""
# from AAComposition import CalculateAAComposition, CalculateDipeptideComposition, GetSpectrumDict

from Autocorrelation import CalculateNormalizedMoreauBrotoAutoTotal, CalculateMoranAutoTotal, CalculateGearyAutoTotal

from Autocorrelation import CalculateEachGearyAuto, CalculateEachMoranAuto, CalculateEachNormalizedMoreauBrotoAuto

# from CTD import CalculateCTD

# from QuasiSequenceOrder import GetSequenceOrderCouplingNumberTotal, GetQuasiSequenceOrder, \
#     GetSequenceOrderCouplingNumberp, GetQuasiSequenceOrderp

# from PseudoAAC import _GetPseudoAAC, GetAPseudoAAC, GetPseudoAAC

# from GetSubSeq import GetSubSequence

# from AAIndex import GetAAIndex1, GetAAIndex23

# from ConjointTriad import CalculateConjointTriad


class PyProtein():
    """
    This GetProDes class aims at collecting all descriptor calcualtion modules into a simple class.

    """
    AALetter = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

    Version = 1.0

    def __init__(self, ProteinSequence=''):
        """
        input a protein sequence
        """
        if len(ProteinSequence) == 0:
            print("You must input a protein sequence when constructing a object. It is a string!")
        else:
            self.ProteinSequence = ProteinSequence

    

    def GetMoreauBrotoAuto(self):
        """
        Normalized Moreau-Broto autocorrelation descriptors (240)

        Usage:

        result = GetMoreauBrotoAuto()
        """
        res = CalculateNormalizedMoreauBrotoAutoTotal(self.ProteinSequence)
        return res

    def GetMoranAuto(self):
        """
        Moran autocorrelation descriptors (240)

        Usage:

        result = GetMoranAuto()
        """
        res = CalculateMoranAutoTotal(self.ProteinSequence)
        return res

    def GetGearyAuto(self):
        """
        Geary autocorrelation descriptors (240)

        Usage:

        result = GetGearyAuto()
        """
        res = CalculateGearyAutoTotal(self.ProteinSequence)
        return res



#####################################################################################################