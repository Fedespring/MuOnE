## For starting:

git init

git add tbAnalyzerNew.C

git commit -m "fist commit"

git push https://github.com/Fedespring/MuOnE.git

enter username

token

e speri che vada



## Calibration:

cd RelativeCalCon/

root -b

.L SkAnalyzerNew.C++ 

SkAnalyzer t("/eos/user/e/elusiani/muone/for_federica/skim_Pilot2025/Ele20GeV_250715.root")

t.iterate(4,-20,100, "Ele20GeV_250715.hist")



## Analysis:

cd Calibrated/

root -b

.L SkAnalyzerNew.C++ 

SkAnalyzer t("/eos/user/e/elusiani/muone/for_federica/skim_Pilot2025/Ele20GeV_250715.root")

t.Loop(-20, 202507, 0, "Ele20GeV_250715.root")