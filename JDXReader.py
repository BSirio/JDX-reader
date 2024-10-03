import matplotlib.pyplot as plt
import numpy as np
import cv2
import os
from scipy import signal
from datetime import datetime, timezone

#filePath  = "C:\\Users\\BRUNIASA\\Desktop\\KAUST\\PhD\\WORKS\\co2 capture\\breakthrough ftir\\75C_50sscm_desorption_plus_dead_volume\\Lewatiti_75C50_sscm_ads_des_jdx.JDX"
#outFolder = "C:\\Users\\BRUNIASA\\Desktop\\KAUST\\PhD\\WORKS\\co2 capture\\breakthrough ftir\\75C_50sscm_desorption_plus_dead_volume"

#filePath  = "C:\\Users\\BRUNIASA\\Desktop\\KAUST\\PhD\\WORKS\\co2 capture\\breakthrough ftir\\25C_50sscm_wet_15C\\Lewatiti_25C50_sscm__ads_des_wet_T_chiller_15C_all.JDX"
#outFolder = "C:\\Users\\BRUNIASA\\Desktop\\KAUST\\PhD\\WORKS\\co2 capture\\breakthrough ftir\\25C_50sscm_wet_15C"

#filePath  = "C:\\Users\\BRUNIASA\\Desktop\\KAUST\\PhD\\WORKS\\co2 capture\\breakthrough ftir\\dead_volume\\dead_volume_wet_40sscm.JDX"
#outFolder = "C:\\Users\\BRUNIASA\\Desktop\\KAUST\\PhD\\WORKS\\co2 capture\\breakthrough ftir\\dead_volume"
#filePath  = "C:\\Users\\BRUNIASA\\Desktop\\KAUST\\PhD\\WORKS\\co2 capture\\breakthrough ftir\\ads-des_25C_Tsorb_15C_Tchiller_40sscm\\ads-des_25C_Tsorb_15C_Tchiller_40sscm.JDX"
#outFolder = "C:\\Users\\BRUNIASA\\Desktop\\KAUST\\PhD\\WORKS\\co2 capture\\breakthrough ftir\\ads-des_25C_Tsorb_15C_Tchiller_40sscm"

filePath  = "C:\\Users\\BRUNIASA\\Desktop\\KAUST\\PhD\\WORKS\\co2 capture\\breakthrough ftir\\zulmas_experiments\\v2\\Purolite__humid_138_mg_ads_25_C_des_25C_RH_60_3_jun_v2.JDX"
outFolder = "C:\\Users\\BRUNIASA\\Desktop\\KAUST\\PhD\\WORKS\\co2 capture\\breakthrough ftir\\zulmas_experiments\\v2"


inFile = open(filePath, 'r')

plotSpectra = False
generateVideo = False

#A = 1/153.0899  # measured values
A = 0.006  # measured values
CO2WaveLims = [2270, 2390]
H2OWaveLims = [1502, 1513]

#find number of plots
numOfPlots = 0
while True:
    line = inFile.readline()
    #print(line[0:8])
    if not line:
        print("ERROR: \"BLOCKS\" keyword not found!")
        exit()
    if line[0:8] == "##BLOCKS":
        words = line.split('=')
        numOfPlots = int(words[1])
        break

print("%s spectra are present in the file" % numOfPlots)

print("Extracting data ...")
# extract the spectra data
spectra = [[[],[]] for i in range(numOfPlots)]
background = [[],[]]
times = [0 for i in range(numOfPlots)]
dateformat = "%a %b %d %H:%M:%S %Y (%Z%z)"

#print(spectra)
index = 0
while index < numOfPlots:
    deltaX = 0
    xFactor = 0
    yFactor = 0
    startDataAquisition = False
    titleFound = False
    backgroundFrame = False
    Xs = []
    Ys = []
    #print(index)
    while True:
        line = inFile.readline()
        if startDataAquisition:
            if line[0:5] == "##END":
                break
            words = line.split()
            startX = float(words[0])
            numY = 0
            for YString in words[1:]:
                Y = float(YString)*yFactor
                X = startX*xFactor + deltaX*numY
                Xs.append(X)
                Ys.append(Y)
                numY = numY+1
            continue

        if not line:
            print("ERROR: did not find all the stated spectra before end of file!")
            exit()
        
        if line[0:7] == "##TITLE":
            words = line.split('=')
            #print(".%s." % words[1])
            if "Back" in words[1]: #words[1] == "Back\n":
                print("Background measurement found")
                numOfPlots = numOfPlots-1
                #index = index - 1
                spectra = spectra[:-1]
                times = times[:-1]
                backgroundFrame = True
            if words[1][0] == '*':
                print("Asterisk measurement found")
                numOfPlots = numOfPlots-1
                index = index - 1
                spectra = spectra[:-1]
                times = times[:-1]
                break
            if not backgroundFrame:
                try:
                    date_obj = datetime.strptime(words[1][:-1], dateformat)
                    times[index] = (date_obj - datetime(1970,1,1,tzinfo=timezone.utc)).total_seconds()
                    #if(times[index] == 0):
                    #    print(words[1])
                    titleFound = True
                except:
                    print("Error %s not convertible" % words[1])
                    continue
        #if line[0:6] == "##TIME":
        #    words = line.split("=")
        #    timeString = words[1]
        #    HMS = timeString.split(':')
        #    time = int(HMS[0])*3600 + int(HMS[1])*60 + int(HMS[2])
        #    times[index] = time
        if(line[0:9] == "##XFACTOR"):
            words = line.split("=")
            xFactor = float(words[1])
        if(line[0:9] == "##YFACTOR"):
            words = line.split("=")
            yFactor = float(words[1])
        if(line[0:8] == "##DELTAX"):
            words = line.split("=")
            deltaX = float(words[1])
        if(line[0:8] == "##XYDATA") and titleFound:
            startDataAquisition = True
        if(line[0:8] == "##XYDATA") and backgroundFrame:
            startDataAquisition = True
    if backgroundFrame:
        background[0] = Xs
        background[1] = Ys
        backgroundFrame = False
        continue
    spectra[index][0] = Xs
    spectra[index][1] = Ys    
    index = index + 1
    #print(index)
inFile.close()
print("    ... COMPLETE!")


print("Processing data ...")
#process the data
minTime = min(times)
#print(minTime)
for index in range(len(times)):
    #print(times[index])
    times[index] = times[index] - minTime
print("Total time = %d s" % (max(times)-min(times)))
# find maximum CO2 absorbance
maxY = 0
for spectrum in spectra:
    for y in spectrum[1]:
        if y > maxY:
            maxY = y



CO2_values = [0 for i in range(numOfPlots)]
for index in range(numOfPlots):
    t = spectra[index][0]
    absorb = spectra[index][1]
    value = 0
    for ind2 in range(len(t)):
        if t[ind2] > CO2WaveLims[0] and t[ind2] < CO2WaveLims[1] and ind2>1:
            value = value + (absorb[ind2] + absorb[ind2-1])*(t[ind2]-t[ind2-1])/2
    CO2_values[index] = value

H2O_values = [0 for i in range(numOfPlots)]
for index in range(numOfPlots):
    t = spectra[index][0]
    absorb = spectra[index][1]
    value = 0
    for ind2 in range(len(t)):
        if t[ind2] > H2OWaveLims[0] and t[ind2] < H2OWaveLims[1] and ind2>1:
            value = value + (absorb[ind2] + absorb[ind2-1])*(t[ind2]-t[ind2-1])/2
    H2O_values[index] = value

# normalize CO2 values
maxCO2value = max(CO2_values)
for index in range(numOfPlots):
    CO2_values[index] = CO2_values[index]/maxCO2value

# normalize H2O values
maxH2Ovalue = max(H2O_values)
for index in range(numOfPlots):
    H2O_values[index] = H2O_values[index]/maxH2Ovalue

# sort vector
times , CO2_values, H2O_values, spectra = (list(t) for t in zip(*sorted(zip(times, CO2_values, H2O_values, spectra))))

#interpolate the concentration profiles
timeStep = 20
numPoints = int(max(times)/timeStep)+1
timesInterp = np.linspace(0, timeStep*(numPoints-1), num = numPoints)
CO2Interp = np.interp(timesInterp, times, CO2_values)
H2OInterp = np.interp(timesInterp, times, H2O_values)
#print(timesInterp)
# compute the system convolution

#A = 1/230  # test values
convolution = [A*np.exp(-A*i*timeStep) for i in range(60)]
#convolution = [1-np.exp(-A*(i+1)*timeStep) for i in range(60)]
CO2Deconv, remainder = signal.deconvolve(CO2Interp,convolution)
H2ODeconv, remainder = signal.deconvolve(H2OInterp,convolution)
deconvTimes = [i*timeStep for i in range(len(CO2Deconv))]
#print(deconvoluted)

print("   ... COMPLETE!")

print("Plotting data ...")

plt.figure()
plt.plot(timesInterp, CO2Interp)
plt.xlabel("Time [s]")
plt.ylabel("CO2 [-]")
plt.title("CO2")
#plt.xlim((0, 10000))
plt.savefig("%s\\CO2.png" % outFolder, dpi=500)
plt.close()

plt.figure()
plt.plot(timesInterp, H2OInterp)
plt.xlabel("Time [s]")
plt.ylabel("H2O [-]")
plt.title("H2O")
#plt.xlim((0, 10000))
plt.savefig("%s\\H2O.png" % outFolder, dpi=500)
plt.close()

plt.figure()
plt.plot(deconvTimes, CO2Deconv)
plt.xlabel("Time [s]")
plt.ylabel("CO2 [-]")
plt.title("CO2")
#plt.xlim((0, 10000))
#plt.ylim((-1, 20))
plt.savefig("%s\\CO2_deconv.png" % outFolder, dpi=500)
plt.close()

plt.figure()
plt.plot(deconvTimes, H2ODeconv)
plt.xlabel("Time [s]")
plt.ylabel("H2O [-]")
plt.title("H2O")
#plt.xlim((0, 10000))
#plt.ylim((-1, 20))
plt.savefig("%s\\H2O_deconv.png" % outFolder, dpi=500)
plt.close()

if len(background[0]) > 0:
    with open("%s\\background.csv" % outFolder, 'w') as bkgFile:
        for i in range(len(background[0])):
            bkgFile.write("%f,%f\n" % (background[0][i], background[1][i]))

with open("%s\\CO2.csv" % outFolder, 'w') as outFile:
    outFile.write("Time [s], CO2[-]\n")
    for i in range(numPoints):
        outFile.write("%f,%f\n" % (timesInterp[i], CO2Interp[i]))

with open("%s\\H2O.csv" % outFolder, 'w') as outFile:
    outFile.write("Time [s], H2O[-]\n")
    for i in range(numPoints):
        outFile.write("%f,%f\n" % (timesInterp[i], H2OInterp[i]))


with open("%s\\CO2_deconvoluted.csv" % outFolder, 'w') as outFile:
    outFile.write("Time [s], CO2[-]\n")
    for i in range(len(CO2Deconv)):
        outFile.write("%f,%f\n" % (deconvTimes[i], CO2Deconv[i]))

with open("%s\\H2O_deconvoluted.csv" % outFolder, 'w') as outFile:
    outFile.write("Time [s], H2O[-]\n")
    for i in range(len(H2ODeconv)):
        outFile.write("%f,%f\n" % (deconvTimes[i], H2ODeconv[i]))


if plotSpectra or generateVideo:
    framesFolderPath = "%s\\frames" % outFolder
    if not os.path.exists(framesFolderPath):
        os.makedirs(framesFolderPath)

    for index in range(numOfPlots):
        outFile = open("%s\\%010d_s.csv" % (framesFolderPath, times[index]), 'w')
        outFile.write("wavelenght[1/cm],absorbance[-]\n")
        for i in range(len(spectra[index][0])):
            outFile.write("%f,%f\n" % (spectra[index][0][i],spectra[index][1][i]))
        outFile.close()
        plt.figure()
        plt.plot(spectra[index][0], spectra[index][1])
        plt.ylim((0, maxY))
        plt.xlabel("Wavelength [1/cm]")
        plt.ylabel("Absorbance [-]")
        plt.title("Time = %d s" % times[index])
        plt.savefig("%s\\%010d_s.png" % (framesFolderPath, times[index]))
        plt.close()
    print("   ... COMPLETE!")

    if generateVideo:
        print("Generate video ...")
        video_name = "%s\\video.avi" % outFolder
        images = [img for img in os.listdir(framesFolderPath) if img.endswith(".png")]
        frame = cv2.imread(os.path.join(framesFolderPath, images[0]))
        height, width, layers = frame.shape
        #fourcc = cv2.VideoWriter.fourcc('M', 'S', 'V', 'C');
        fourcc = 0
        video = cv2.VideoWriter(video_name, fourcc, 60, (width,height))

        for image in images:
            video.write(cv2.imread(os.path.join(framesFolderPath, image)))

        cv2.destroyAllWindows()
        video.release()

print("   ... COMPLETE!")