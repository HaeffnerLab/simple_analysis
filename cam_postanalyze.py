from cct.scripts.experiments.Camera.ion_state_detector_1d import ion_state_detector
import numpy as np
import pylab as pl

def get_statelist(timestr, md, rd):
    """Performs post analysis
    timestr: Time string
    md: MeasDay object
    rd: ReadData object
    Reanalyzes camera data and performs a fit for ion position at every data point.
    works only for a single qubit at the moment
    """
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

    imag_list = rd.get_images(timestr)

    pos = []
    statelist = []
    for image in imag_list:
        mean_image = np.mean(image,0)
        fitter = ion_state_detector(positions)
        result, params = fitter.guess_parameters_and_fit(xx, yy, mean_image)
        pos.append(params['pos0'].value)
        fitter.params['amplitude'].value = ion_params['fit_amplitude']
        state = fitter.state_detection(image[:,0,:])[0]
        statelist.append(1-np.mean(state[:,0]))
    return statelist
