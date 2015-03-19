from simple_analysis import get_data as gd
from cct.scripts.experiments.Camera.ion_state_detector_1d import ion_state_detector
import numpy as np

md = gd.MeasDay('2015Mar17')
rd = gd.ReadData('2015Mar17', experiment='RamseyScanGapContrast')

timestr = '2028_20'
params = md.param_dict[timestr]
ion_params = params.IonsOnCamera
ion_params['ion_positions'] = ion_params['ion_positions'][1]
ion_params.vertical_max = ion_params.vertical_min
positions = ion_params['ion_positions']


image_region = [int(ion_params.horizontal_bin),
                int(ion_params.vertical_bin),
                int(ion_params.horizontal_min),
                int(ion_params.horizontal_max),
                int(ion_params.vertical_min),
                int(ion_params.vertical_max)]

x_axis = np.arange(ion_params.horizontal_min, ion_params.horizontal_max + 1, image_region[0])
y_axis = np.arange(ion_params.vertical_min, ion_params.vertical_max + 1, image_region[1])
xx, yy = np.meshgrid(x_axis, y_axis)

fh = open('images.npy')
n_imag = 0
imag_list = []
try:
    while(True):
        a = np.load(fh)
        imag_list.append(a)
except:
    pass
pos = []
for image in imag_list:
    image = image[0,:,:]
    fitter = ion_state_detector(positions)
    result, params = fitter.guess_parameters_and_fit(xx, yy, image)
    pos.append(params['pos0'].value)
