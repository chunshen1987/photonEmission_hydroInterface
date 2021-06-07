#!/usr/bin/env python3

from sys import argv
from os import path
import matplotlib.pyplot as plt
from numpy import *

ViscousList = ["eq", "vis", "bulkvis", "tot",
               "vis_deltaf_restricted", "bulkvis_deltaf_restricted"]

HGList = ["HG_2to2_meson_total", "HG_omega", "HG_pipi_bremsstrahlung",
          "HG_rho_spectralfun"]

QGPList = ["QGP_2to2_total", "QGP_AMYcollinear"]

try:
    RESULTS = path.abspath(argv[1])
except:
    print("Usage: {} RESULTFolder".format(argv[0]))
    exit(0)

visType = 3
print("Making plots from {} ...".format(RESULTS))

npT = 0
# plot photon spectra
fig = plt.figure()
ax = plt.axes([0.12, 0.12, 0.83, 0.83])
for ichannel in QGPList + HGList:
    data = loadtxt(path.join(
        RESULTS, "{}_SPvn_{}.dat".format(ichannel, ViscousList[visType])))
    npT = len(data[:, 0])
    plt.plot(data[:, 0], data[:, 1], label=ichannel)
plt.legend(loc=0)
plt.yscale("log")
plt.xlabel(r"$p_T$ (GeV)")
plt.ylabel(r"$dN/(2\pi p_T dp_T dy)$ (GeV$^{-2}$)")
plt.savefig("photon_spectra_individual_channels.pdf")

pTArr = zeros(npT)
QGPtotal = zeros(npT)
HGtotal = zeros(npT)
Thermaltotal = zeros(npT)
for ichannel in QGPList:
    data = loadtxt(path.join(
        RESULTS, "{}_SPvn_{}.dat".format(ichannel, ViscousList[visType])))
    pTArr = data[:, 0]
    QGPtotal += data[:, 1]
for ichannel in HGList:
    data = loadtxt(path.join(
        RESULTS, "{}_SPvn_{}.dat".format(ichannel, ViscousList[visType])))
    HGtotal += data[:, 1]
Thermaltotal = QGPtotal + HGtotal
fig = plt.figure()
ax = plt.axes([0.12, 0.12, 0.83, 0.83])
plt.plot(pTArr, QGPtotal, label="QGP")
plt.plot(pTArr, HGtotal, label="HG")
plt.plot(pTArr, Thermaltotal, label="Thermal")
plt.legend(loc=0)
plt.yscale("log")
plt.xlabel(r"$p_T$ (GeV)")
plt.ylabel(r"$dN/(2\pi p_T dp_T dy)$ (GeV$^{-2}$)")
plt.savefig("photon_spectra_QGP_and_HG.pdf")

# plot photon v2
fig = plt.figure()
ax = plt.axes([0.12, 0.12, 0.83, 0.83])
for ichannel in QGPList + HGList:
    data = loadtxt(path.join(
        RESULTS, "{}_SPvn_{}.dat".format(ichannel, ViscousList[visType])))
    plt.plot(data[:, 0], data[:, 7], label=ichannel)
plt.legend(loc=0)
plt.xlabel(r"$p_T$ (GeV)")
plt.ylabel(r"$v_2$")
plt.savefig("photon_v2_individual_channels.pdf")

pTArr = zeros(npT)
QGPv2 = zeros(npT)
HGv2 = zeros(npT)
Thermalv2 = zeros(npT)
fig = plt.figure()
ax = plt.axes([0.12, 0.12, 0.83, 0.83])
for ichannel in QGPList:
    data = loadtxt(path.join(
        RESULTS, "{}_SPvn_{}.dat".format(ichannel, ViscousList[visType])))
    pTArr = data[:, 0]
    QGPv2 += data[:, 1]*data[:, 7]
QGPv2 = QGPv2/QGPtotal
for ichannel in HGList:
    data = loadtxt(path.join(
        RESULTS, "{}_SPvn_{}.dat".format(ichannel, ViscousList[visType])))
    HGv2 += data[:, 1]*data[:, 7]
HGv2 = HGv2/HGtotal
Thermalv2 = (QGPv2*QGPtotal + HGv2*HGtotal)/Thermaltotal
fig = plt.figure()
ax = plt.axes([0.12, 0.12, 0.83, 0.83])
plt.plot(pTArr, QGPv2, label="QGP")
plt.plot(pTArr, HGv2, label="HG")
plt.plot(pTArr, Thermalv2, label="Thermal")
plt.legend(loc=0)
plt.xlabel(r"$p_T$ (GeV)")
plt.ylabel(r"$v_2$")
plt.savefig("photon_v2_QGP_vs_HG.pdf")
