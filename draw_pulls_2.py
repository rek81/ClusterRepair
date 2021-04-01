import math
import ROOT
from ROOT import *
import os


F = TFile("../../python/PixelTree_DamageclusterFalse.root")
T = F.Get("pixelTree")
n = T.GetEntries()


for layer in range(1,4):
        

        px = TH1F("h_x_L"+str(layer), ";"+str(layer)+";Events", 100, -300, 300)
        py = TH1F("h_y_L"+str(layer), ";"+str(layer)+";Events", 100, -300, 300)
        ppullx = TH1F("h_xpull_L"+str(layer), ";"+str(layer)+";Events", 100, -3, 3)
        ppully = TH1F("h_ypull_L"+str(layer), ";"+str(layer)+";Events", 100, -3, 3)
        
        for i in range(n):
#                print i
                
                T.GetEntry(i)
                
                if T.ClN<20000:

                        print T.ClN
                        for j in T.ClN:
#                                print j
                                
                                if T.ClType[j]==1: 
                                        
                                        for k in T.ClSimHitN[j]:

                                                if T.ClLayer[j]==layer and T.ClRhHasBadPixels[j]>0.:
#                                                if T.ClLayer[j]==layer and T.ClRhHasBadPixels[j]>0. and abs(T.ClSimHitLx[k]-T.ClRhLx[j])<9999.:
                                                        
                                                        xres = 100.*(T.ClRhLx[j] - T.ClSimHitLx[k])
                                                        yres = 100.*(T.ClRhLy[j] - T.ClSimHitLy[k])
                                                        xpull = ((T.ClRhLx[j] - T.ClSimHitLx[k])/T.ClRhLxE[j])
                                                        ypull = ((T.ClRhLy[j] - T.ClSimHitLy[k])/T.ClRhLyE[j])
                                                        
                                                        px.Fill(xres)
                                                        py.Fill(yres)
                                                        ppullx.Fill(xpull)
                                                        ppully.Fill(ypull)
        

        c1b = TCanvas("c1_L"+str(layer), "", 0, 0, 800, 800)
        px.Draw("pe")
        px.Fit("gaus")
#        px.GetFunction("gaus").SetLineColor(kBlue)
        c1b.Print("badpix_residx_False_L"+str(layer)+".png")
        
        c2b = TCanvas("c2_L"+str(layer), "", 0,0 , 800, 800)
        py.Draw("pe")
        py.Fit("gaus")
#        py.GetFunction("gaus").SetLineColor(kBlue)
        c2b.Print("badpix_residy_False_L"+str(layer)+".png")
        
        c3b = TCanvas("c3_L"+str(layer), "", 0, 0, 800, 800)
        ppullx.Draw("pe")
        ppullx.Fit("gaus")
#        ppullx.GetFunction("gaus").SetLineColor(kBlue)
        c3b.Print("badpix_pullx_False_L"+str(layer)+".png")
        
        c4b = TCanvas("c4_L"+str(layer), "", 0,0 , 800, 800)
        ppully.Draw("pe")
        ppully.Fit("gaus")
#        ppully.GetFunction("gaus").SetLineColor(kBlue)
        c4b.Print("badpix_pully_False_L"+str(layer)+".png")

