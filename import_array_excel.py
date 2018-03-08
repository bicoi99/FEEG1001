
from aerofoil_optimization import parametric_aerofoil_four_var
import xlsxwriter as xlsx
import numpy as np

w = np.array([0.51223816,  0.64842888,  0.61496199,  0.47028275])
points = parametric_aerofoil_four_var(w)

wb = xlsx.Workbook('spline_data.xlsx')
ws = wb.add_worksheet()
row = 0
col = 0

for x, y in points:
    ws.write(row, col, x * 100)
    ws.write(row, col + 1, y * 100)
    ws.write(row, col + 2, 0)
    row += 1

wb.close()
