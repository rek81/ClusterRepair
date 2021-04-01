import ROOT
import math
import os

RDF = ROOT.ROOT.RDataFrame
ROOT.ROOT.EnableImplicitMT() 

#auto cotGcut = [](float x){return 1/math.cot(TkBeta)>0.}

for layer in range(1,4):


        F = RDF("pixelTree", "../../python/PixelTree_DamageclusterFalse.root")
        rdf = F.Define("total_weight", "1.0")

        h = rdf.Filter("ClType==1 && ClN<20000 && "+"ClLayer=="+str(layer)+" && ClRhHasBadPixels>0")
 #       hcotG = rdf.Filter("&& ClType==1. && ClN<20000&&"+"ClLayer=="+layer+". && ClRhHasBadPix>0. && 1/TMath.Tan()")
 #       hcotL = rdf.Filter("&& ClType==1. && ClN<20000&&"+"ClLayer=="+layer+". && ClRhHasBadPix>0.")

#        p_x = h.Histo1D(("h_x_L"+str(layer), ";"+str(layer)+";Events", 100, -300, 300), "ClRhLx-ClSimHitLx", "total_weight")
        p_x = h.Histo1D(("h_x_L"+str(layer), ";"+str(layer)+";Events", 100, -300, 300), "ClRhLx", "total_weight")
        px = p_x.DrawCopy()
        pxdraw = px.Clone("px")
        p_y = h.Histo1D(("h_y_L"+str(layer), ";"+str(layer)+";Events", 100, -300, 300), "ClRhLy-ClSimHitLy", "total_weight")
        py = p_y.DrawCopy()
        pydraw = py.Clone("py")
        p_pullx = h.Histo1D(("h_xpull_L"+str(layer), ";"+str(layer)+";Events", 100, -3, 3), "(ClRhLx-ClSimHitLx)/ClRhLxE", "total_weight")
        ppullx = p_pullx.DrawCopy()
        ppullxdraw = px.Clone("ppullx")
        p_pully = h.Histo1D(("h_ypull_L"+str(layer), ";"+str(layer)+";Events", 100, -3, 3), "(ClRhLy-ClSimHitLy)/ClRhLyE", "total_weight")
        ppully = p_pully.DrawCopy()
        ppullydraw = ppully.Clone("ppully")
        

        c1b = TCanvas("c1_L"+str(layer), "", 0, 0, 800, 800)
        pxdraw.Draw("pe")
        pxdraw.Fit("gaus")
        pxdraw.GetFunction("gaus").SetLineColor(kBlue);
        c1b.Print("badpix_residx_False_L"+str(layer)+".png")
        
        c2b = TCanvas("c2_L"+str(layer), "", 0,0 , 800, 800)
        pydraw.Draw("pe")
        pydraw.Fit("gaus")
        pydraw.GetFunction("gaus").SetLineColor(kBlue)
        c2b.Print("badpix_residy_False_L"+str(layer)+".png")
        
        c3b = TCanvas("c3_L"+str(layer), "", 0, 0, 800, 800)
        ppullxdraw.Draw("pe")
        ppullxdraw.Fit("gaus")
        ppullxdraw.GetFunction("gaus").SetLineColor(kBlue)
        c3b.Print("badpix_pullx_False_L"+str(layer)+".png")
        
        c4b = TCanvas("c4_L"+str(layer), "", 0,0 , 800, 800)
        ppullydraw.Draw("pe")
        ppullydraw.Fit("gaus")
        ppullydraw.GetFunction("gaus").SetLineColor(kBlue)
        c4b.Print("badpix_pully_False_L"+str(layer)+".png")

