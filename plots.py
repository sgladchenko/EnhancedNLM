#!/usr/bin/env python3

import numpy as np
import re

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import ticker
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

from scipy.fftpack import fft, ifft, fftfreq
from scipy.optimize import curve_fit as scfit

import numpy as np

tau1 = np.array([[0 + 0j, 1 + 0j], [1 + 0j, 0 + 0j]])
tau2 = np.array([[0 + 0j, 0 - 1j], [0 + 1j, 0 + 0j]])
tau3 = np.array([[1 + 0j, 0 + 0j], [0 + 0j, -1 + 0j]])

def Fix_E(filename, filename_zgrid, filename_xgrid, filename_plot, i=0, j=0):
	data  = np.fromfile(filename, dtype=np.complex128)
	zgrid = np.fromfile(filename_zgrid, dtype=np.float64)
	xgrid = np.fromfile(filename_xgrid, dtype=np.float64)

	N_Z = len(zgrid)
	N_X = len(xgrid)

	#mesh = np.array([[np.real(data[z*N_X*4 + x*4 + i*2 + j]) for x in range(N_X)] + [np.real(data[z*N_X*4 + x*4 + i*2 + j]) for x in range(N_X)] for z in range(N_Z)])
	mesh = np.array([[np.real(data[z*N_X*4 + x*4 + i*2 + j]) for x in range(N_X)] for z in range(N_Z)])

	Zs = np.array(zgrid)
	#Xs = np.array(list(xgrid) + [x + xgrid[-1] for x in xgrid])
	Xs = np.array(xgrid)

	fig = Figure(figsize=(8, 6))
	FigureCanvas(fig)

	Xs, Zs = np.meshgrid(Xs, Zs)

	axs  = fig.add_subplot(111)
	plot = axs.pcolor(Xs, Zs, mesh, cmap=cm.viridis, vmin=0, vmax=1)

	axs.set_title(r"$\rho$[{}, {}]".format(i, j), fontsize=12)
	axs.set_xlabel(r"$x, \mathrm{km}$")
	axs.set_ylabel(r"$z, \mathrm{km}$")

	fig.colorbar(plot)
	fig.savefig(filename_plot, fmt="png")

def Average_by_x(filename, filename_zgrid, filename_egrid, filename_plot, i=0, j=0):
	data  = np.fromfile(filename, dtype=np.complex128)
	zgrid = np.fromfile(filename_zgrid, dtype=np.float64)
	egrid = np.fromfile(filename_egrid, dtype=np.float64)

	N_Z = len(zgrid)
	N_E = len(egrid)

	mesh = np.array([[np.real(data[z*N_E*4 + e*4 + i*2 + j]) for e in range(N_E)] for z in range(N_Z)])

	Zs = np.array(zgrid)
	Es = np.array(egrid)

	fig = Figure(figsize=(8, 6))
	FigureCanvas(fig)

	Es, Zs = np.meshgrid(Es, Zs)

	axs  = fig.add_subplot(111)
	plot = axs.pcolor(Es, Zs, mesh, cmap=cm.viridis, vmin=0, vmax=1)

	axs.set_title(r"$\rho$[{}, {}]".format(i, j), fontsize=12)
	axs.set_xlabel(r"$E, \mathrm{MeV}$")
	axs.set_ylabel(r"$z, \mathrm{km}$")

	fig.colorbar(plot)
	fig.savefig(filename_plot, fmt="png")

def Average_by_x_curves(filename, filename_zgrid, filename_egrid, filename_plot, i=0, j=0):
	data  = np.fromfile(filename, dtype=np.complex128)
	zgrid = np.fromfile(filename_zgrid, dtype=np.float64)
	egrid = np.fromfile(filename_egrid, dtype=np.float64)

	N_Z = len(zgrid)
	N_E = len(egrid)

	mesh = np.array([[np.real(data[z*N_E*4 + e*4 + i*2 + j]) for z in range(N_Z)] for e in range(N_E)])

	fig = Figure(figsize=(8, 6))
	FigureCanvas(fig)

	axs  = fig.add_subplot(111)
	fig.subplots_adjust(right=0.80)

	errorbars = [None for i in range(N_E)]
	colors   = [cm.viridis(1.0 - e/N_E) for e in range(N_E)]

	for e in range(N_E):
		errorbars[e] = axs.errorbar(
			x=zgrid,
			y=mesh[e],
			fmt='-',
			ecolor=colors[e],
			color=colors[e],
			markersize=6
		)			

	axs.legend(errorbars, ["E = {:.2f} MeV".format(1.0 + e*49/N_E) for e in range(N_E)], fontsize=6, loc="lower right", bbox_to_anchor=(1.25, 0.0))

	axs.set_title(r"$\rho$[{}, {}]".format(i, j), fontsize=12)
	axs.set_xlabel(r"$z, \mathrm{km}$")
	axs.set_xlim([np.min(zgrid), np.max(zgrid)])
	axs.set_ylim([0, 1])

	fig.savefig(filename_plot, fmt="png")

def Fourier_curves(filename_rec, filename_xgrid, filename_egrid, filename_plot, func=0, i=0, j=0):
	data  = np.fromfile(filename_rec, dtype=np.complex128)
	xgrid = np.fromfile(filename_xgrid, dtype=np.float64)
	egrid = np.fromfile(filename_egrid, dtype=np.float64)

	N_X = len(xgrid)
	N_E = len(egrid)

	arrs = np.array([[np.real(data[x*N_E*16 + e*16 + func*4 + i*2 + j]) for x in range(N_X)] for e in range(N_E)])

	#arrs = [[np.sin(2*np.pi*n/N_X) for n in range(N_X)], ]

	fig = Figure(figsize=(8, 6))
	FigureCanvas(fig)

	axs  = fig.add_subplot(111)
	fig.subplots_adjust(right=0.80)

	errorbars = [None for i in range(N_E)]
	colors   = [cm.viridis(1.0 - e/N_E) for e in range(N_E)]

	Fs = None

	for e in range(N_E):#(5, 15, 25):
		As = fft(arrs[e])
		Fs = fftfreq(N_X, d=1/N_X)

		errorbars[e] = axs.errorbar(
			x=[x for x in Fs if x > 0],
			y=[np.abs(As[n]) for n in range(N_X) if Fs[n] > 0],
			fmt='-',
			ecolor=colors[e],
			color=colors[e],
			markersize=6
		)			

	axs.legend(errorbars, ["E = {:.2f} MeV".format(1.0 + e*49/N_E) for e in range(N_E)], fontsize=6, loc="lower right", bbox_to_anchor=(1.25, 0.0))

	axs.set_title(r"$\rho$[{}, {}]".format(i, j), fontsize=12)
	axs.set_xlabel(r"$k$")
	axs.set_ylabel(r"$\vert C_k \vert$")
	#axs.set_xlim([0, 400])
	axs.set_yscale('log')

	fig.savefig(filename_plot, fmt="png")

def Last_z_curves(filename_rec, filename_xgrid, filename_egrid, filename_plot, func=0, i=0, j=0):
	data  = np.fromfile(filename_rec, dtype=np.complex128)
	xgrid = np.fromfile(filename_xgrid, dtype=np.float64)
	egrid = np.fromfile(filename_egrid, dtype=np.float64)

	N_X = len(xgrid)
	N_E = len(egrid)

	arrs = np.array([[np.real(data[x*N_E*16 + e*16 + func*4 + i*2 + j]) for x in range(N_X)] for e in range(N_E)])

	fig = Figure(figsize=(8, 6))
	FigureCanvas(fig)

	axs  = fig.add_subplot(111)
	fig.subplots_adjust(right=0.80)

	errorbars = [None for i in range(N_E)]
	colors   = [cm.viridis(1.0 - e/N_E) for e in range(N_E)]

	Fs = None

	for e in (5, 15, 25):
		errorbars[e] = axs.errorbar(
			x=xgrid,
			y=arrs[e],
			fmt='-',
			ecolor=colors[e],
			color=colors[e],
			markersize=6
		)			

	#axs.legend(errorbars, ["E = {:.2f} MeV".format(1.0 + e*49/N_E) for e in (5, 15, 25)], fontsize=6, loc="lower right", bbox_to_anchor=(1.25, 0.0))

	axs.set_title(r"$\rho$[{}, {}]".format(i, j), fontsize=12)
	axs.set_xlabel(r"$x, \mathrm{km}$")
	axs.set_xlim([np.min(xgrid), np.max(xgrid)])

	fig.savefig(filename_plot, fmt="png")

def Fourier_curves_Pauli_components(filename_rec, filename_xgrid, filename_egrid, filename_plot, func=0, max_k=350, loglog=True):
	data  = np.fromfile(filename_rec, dtype=np.complex128)
	xgrid = np.fromfile(filename_xgrid, dtype=np.float64)
	egrid = np.fromfile(filename_egrid, dtype=np.float64)

	N_X = len(xgrid)
	N_E = len(egrid)

	mem = [[np.array([[0.0 + 0.0j, 0.0 + 0.0j], [0.0 + 0.0j, 0.0 + 0.0j]]) for x in range(N_X)] for e in range(N_E)]

	for e in range(N_E):
		for x in range(N_X):
			for i in range(2):
				for j in range(2):
					mem[e][x][i][j] = data[x*N_E*16 + e*16 + func*4 + i*2 + j]

	data = [[[0.0 + 0.0j for x in range(N_X)] for e in range(N_E)] for index in range(3)]

	for index, tau in enumerate((tau1, tau2, tau3)):
		for e in range(N_E):
			for x in range(N_X):
				tmp = np.dot(np.matrix(mem[e][x]).getH(), np.matrix(tau))
				data[index][e][x] = 0.5 * np.trace(tmp)

	fig = Figure(figsize=(10, 10))
	FigureCanvas(fig)
	fig.subplots_adjust(right=0.80)

	axs1  = fig.add_subplot(311)
	axs2  = fig.add_subplot(312)
	axs3  = fig.add_subplot(313)

	es = range(N_E)
	errorbars = []

	for index, axs in enumerate((axs1, axs2, axs3)):
		errorbars = [None for i in range(len(es))]
		colors   = [cm.viridis(1.0 - es[i]/N_E) for i in range(len(es))]

		Fs = None

		for counter, e in enumerate(es):
			As = fft(data[index][e])
			Fs = fftfreq(N_X, d=1/N_X)

			errorbars[counter] = axs.errorbar(
				x=[x for x in Fs if x > 0],
				y=[np.abs(As[n]) for n in range(N_X) if Fs[n] > 0],
				fmt='-',
				ecolor=colors[counter],
				color=colors[counter],
				markersize=6
			)			

		axs.set_ylabel(r"$\vert \rho_k \vert$, Component: $\tau_{" + "{}".format(index + 1) + r"}$")
		axs.set_xlim([1, max_k])
		axs.set_ylim([1e-12, 1e4])
		axs.grid(True)
		axs.set_yscale('log')

		if loglog:
			axs.set_xscale('log')

	axs3.set_xlabel(r"$k$")

	axs3.legend(errorbars, ["E = {:.2f} MeV".format(1.0 + es[counter]*49/N_E) for counter in range(len(es))],
			    title="Spectral components", fontsize=7, loc="lower right", bbox_to_anchor=(1.25, 0.0)
	)

	fig.savefig(filename_plot, fmt="png", bbox_inches='tight')

def Last_z_curves_Pauli_components(filename_rec, filename_xgrid, filename_egrid, filename_plot, func=0):
	data  = np.fromfile(filename_rec, dtype=np.complex128)
	xgrid = np.fromfile(filename_xgrid, dtype=np.float64)
	egrid = np.fromfile(filename_egrid, dtype=np.float64)

	N_X = len(xgrid)
	N_E = len(egrid)

	mem = [[np.array([[0.0 + 0.0j, 0.0 + 0.0j], [0.0 + 0.0j, 0.0 + 0.0j]]) for x in range(N_X)] for e in range(N_E)]

	for e in range(N_E):
		for x in range(N_X):
			for i in range(2):
				for j in range(2):
					mem[e][x][i][j] = data[x*N_E*16 + e*16 + func*4 + i*2 + j]

	data = [[[0.0 + 0.0j for x in range(N_X)] for e in range(N_E)] for index in range(3)]

	for index, tau in enumerate((tau1, tau2, tau3)):
		for e in range(N_E):
			for x in range(N_X):
				tmp = np.dot(np.matrix(mem[e][x]).getH(), np.matrix(tau))
				data[index][e][x] = 0.5 * np.trace(tmp)

	fig = Figure(figsize=(10, 10))
	FigureCanvas(fig)
	fig.subplots_adjust(right=0.80)

	axs1  = fig.add_subplot(311)
	axs2  = fig.add_subplot(312)
	axs3  = fig.add_subplot(313)

	es = range(N_E)
	errorbars = []

	for index, axs in enumerate((axs1, axs2, axs3)):
		errorbars = [None for i in range(len(es))]
		colors   = [cm.viridis(1.0 - es[i]/N_E) for i in range(len(es))]

		for counter, e in enumerate(es):
			errorbars[counter] = axs.errorbar(
				x=xgrid,
				y=np.abs(data[index][e]),
				fmt='-',
				ecolor=colors[counter],
				color=colors[counter],
				markersize=6
			)			

		axs.set_ylabel(r"$|\rho_{" + "{}".format(index + 1) + r"}(x)|$, $\zeta=+1$")
		axs.set_xlim([np.min(xgrid), np.max(xgrid)])
		axs.set_ylim([0, 1.0])
		axs.grid(True)

	axs3.legend(errorbars, ["E = {:.2f} MeV".format(1.0 + es[counter]*49/N_E) for counter in range(len(es))],
			    title="Spectral components", fontsize=7, loc="lower right", bbox_to_anchor=(1.25, 0.0)
	)

	axs3.set_xlabel(r"$x, \mathrm{km}$")
	fig.savefig(filename_plot, fmt="png", bbox_inches='tight')

# Show the plot, depicting the adiabaticity at each point

def Ad(filename, filename_zgrid, filename_xgrid, filename_egrid, filename_plot, e = 0):
	data  = np.fromfile(filename, dtype=np.complex128)
	zgrid = np.fromfile(filename_zgrid, dtype=np.float64)
	xgrid = np.fromfile(filename_xgrid, dtype=np.float64)
	egrid = np.fromfile(filename_egrid, dtype=np.float64)

	N_Z = len(zgrid)
	N_X = len(xgrid)
	N_E = len(egrid)

	grid = np.array([[np.real(data[z*N_X*N_E + x*N_E + e]) for x in range(N_X)] for z in range(N_Z)])

	Zs = np.array(zgrid)
	Xs = np.array(xgrid)

	fig = Figure(figsize=(8, 6))
	FigureCanvas(fig)

	Xs, Zs = np.meshgrid(Xs, Zs)

	axs  = fig.add_subplot(111)
	plot = axs.pcolor(Xs, Zs, grid, cmap=cm.viridis, vmin=0, vmax=40)

	axs.set_title(r"Adiabaticity", fontsize=12)
	axs.set_xlabel(r"$x, \mathrm{km}$")
	axs.set_ylabel(r"$z, \mathrm{km}$")

	fig.colorbar(plot)
	fig.savefig(filename_plot, fmt="png")

def Ad_one_line(filename, filename_zgrid, filename_xgrid, filename_egrid, filename_plot, e, x):
	ticker.rcParams['xtick.direction'] = 'in'
	ticker.rcParams['ytick.direction'] = 'in'

	data  = np.fromfile(filename, dtype=np.complex128)
	zgrid = np.fromfile(filename_zgrid, dtype=np.float64)
	xgrid = np.fromfile(filename_xgrid, dtype=np.float64)
	egrid = np.fromfile(filename_egrid, dtype=np.float64)

	N_Z = len(zgrid)
	N_X = len(xgrid)
	N_E = len(egrid)

	grid = np.array([np.real(data[z*N_X*N_E + x*N_E + e]) for z in range(N_Z)])

	fig = Figure(figsize=(8, 6))
	FigureCanvas(fig)

	axs  = fig.add_subplot(111)
	fig.subplots_adjust(right=0.80)

	errorbar = axs.errorbar(
			x=zgrid,
			y=grid,
			fmt='-',
			markersize=6
		)			

	axs.set_title(r"Adiabaticity $x={}$ of $N_X={}$".format(x, N_X), fontsize=12)
	axs.set_xlabel(r"$z, \mathrm{km}$")
	axs.set_xlim([0, np.max(zgrid)])
	axs.set_ylim([0, 10])
	axs.grid(True)

	fig.savefig(filename_plot, fmt="png")

def Ad_integral(filename, filename_zgrid, filename_xgrid, filename_egrid, filename_plot, e, x):
	ticker.rcParams['xtick.direction'] = 'in'
	ticker.rcParams['ytick.direction'] = 'in'

	data  = np.fromfile(filename, dtype=np.complex128)
	zgrid = np.fromfile(filename_zgrid, dtype=np.float64)
	xgrid = np.fromfile(filename_xgrid, dtype=np.float64)
	egrid = np.fromfile(filename_egrid, dtype=np.float64)

	N_Z = len(zgrid)
	N_X = len(xgrid)
	N_E = len(egrid)

	grid0 = np.array([np.real(data[z*N_X*N_E + x*N_E + e]) for z in range(N_Z)])
	grid = []
	tmp = 0

	for z in range(len(zgrid) - 1):
		tmp += 1 / grid0[z] * (zgrid[z+1] - zgrid[z]);
		grid.append(tmp)

	fig = Figure(figsize=(8, 6))
	FigureCanvas(fig)

	axs  = fig.add_subplot(111)
	fig.subplots_adjust(right=0.80)

	errorbar = axs.errorbar(
			x=zgrid[:-1],
			y=grid,
			fmt='-',
			markersize=6
		)			

	axs.set_title(r"Integral adiabaticity $x={}$ of $N_X={}$".format(x, N_X), fontsize=12)
	axs.set_xlabel(r"$z, \mathrm{km}$")
	axs.set_xlim([0, np.max(zgrid)])
	#axs.set_ylim([0, 10])
	axs.grid(True)

	fig.savefig(filename_plot, fmt="png")

if __name__ == "__main__":
	Fourier_curves_Pauli_components("./data/rec.bin", "./data/xgrid.bin", "./data/egrid.bin", "./plots/Fourier_curves_Pauli_components.png", max_k=350, loglog=False)
	Last_z_curves_Pauli_components("./data/rec.bin", "./data/xgrid.bin", "./data/egrid.bin", "./plots/Last_z_curves_Pauli_components.png")

	Fix_E("./data/Fix_E.bin", "./data/zgrid.bin", "./data/xgrid.bin", "./plots/Fix_E.png")
	Average_by_x("./data/Average_by_x.bin", "./data/zgrid.bin", "./data/egrid.bin", "./plots/Average_by_x.png")
	Average_by_x_curves("./data/Average_by_x.bin", "./data/zgrid.bin", "./data/egrid.bin", "./plots/Average_by_x_curves.png")
