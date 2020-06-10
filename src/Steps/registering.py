# -*- coding: utf-8 -*-

"""

Created on Tue Feb  4 14:47:00 2020

@author: Melisa, Morgane

"""

import math
import os
import pickle

import configuration
import numpy as np
from caiman.base.rois import register_multisession
from caiman.source_extraction.cnmf.cnmf import load_CNMF

import Steps.normalized_traces as normalization
from Database.database_connection import database

cursor = database.cursor()


# Method possibilities (model method): registration (True) or matching (False)
# cost_threshold: threshold for cost in matching with Hungarian matching algorithm.
# max_dist : maximum distance between centroids to allow a matching.
# max_cell_size and min_cell size should be taken from the distribution of typical sizes (use function typical size)
# parameters = { 'session_wise': False,'model_method': False, 'cost_threshold' : 0.9 , 'max_dist' : 15 ,
#               'min_cell_size' : 10, 'max_cell_size' : 25}

class estimates(object):
    def __init__(self, A=None, C=None):
        self.A = A
        self.C = C


def run_registration(mouse, session, component_evaluation_v):
    """
    This is the main registering function. Is is supposed to be run after trial wise component evaluation.
    Registration takes over different contours of trial wise source extracted contours and do a matching between cells.
    It can use two different methods: Hungarian matching algorithm (RegisterMulti) (as implement in Giovannucci, et al.
    2019) or cell registration (CellReg)using centroids distance and spatial correlation (as implemented in Sheintuch, et al. 2017).
    Default method is registration with no modeling of distributions of centroids and spatial correlation.

    Args:
        component_evaluation_v: version of the previous step
        mouse: mouse's number
        session: which session you want to analyse
    """

    # Sort the dataframe correctly

    sql = "SELECT component_evaluation_main  FROM Analysis WHERE mouse = ? AND session=? AND component_evaluation_v=?"
    val = [mouse, session, component_evaluation_v]
    cursor.execute(sql, val)
    result = cursor.fetchall()
    file = []
    inter = []
    for x in result:
        inter += x
    for y in inter:
        file.append(y)

    trials_order = []
    for i in range(len(file)):
        sql = "SELECT is_rest FROM Analysis WHERE component_evaluation_main=?"
        val = [file[i]]
        cursor.execute(sql, val)
        result = cursor.fetchall()
        is_rest = []
        trial = []
        for x in result:
            is_rest += x
        sql = "SELECT trial FROM Analysis WHERE component_evaluation_main=?"
        val = [file[i]]
        cursor.execute(sql, val)
        result = cursor.fetchall()
        for x in result:
            trial += x

        for j in range(len(trial)):
            if is_rest[j] == 0:
                order = 2 * trial[j] - 1
            else:
                order = 2 * trial[j]

            sql1 = "UPDATE Analysis SET registration_trials_orders=? WHERE trial=? AND is_rest=? AND mouse = ? AND session=? AND component_evaluation_v=? AND component_evaluation_main=?"
            val1 = [order, trial[j], is_rest[j], mouse, session, component_evaluation_v, file[i]]
            cursor.execute(sql1, val1)
            database.commit()
            trials_order.append(order)

        sql = "SELECT session_wise FROM Analysis WHERE component_evaluation_main=?"
        val = [file[i]]
        cursor.execute(sql, val)
        result = cursor.fetchall()
        session_wise = []
        for x in result:
            session_wise += x
        if session_wise[1] == 'False':
            data_dir = os.environ['DATA_DIR'] + 'data/interim/registration/trial_wise/'
        else:
            data_dir = os.environ['DATA_DIR'] + 'data/interim/registration/session_wise/'

        # determine the file name
        sql = "SELECT mouse,session,trial,is_rest,decoding_v,cropping_v,motion_correction_v,alignment_v,equalization_v,source_extraction_v, component_evaluation_v,registration_v,input,home_path,decoding_main FROM Analysis WHERE component_evaluation_main=?"
        val = [file[i]]
        cursor.execute(sql, val)
        result = cursor.fetchall()
        data = []
        inter = []
        for x in result:
            inter = x
        for y in inter:
            data.append(y)

        # Update the database

        if data[11] == 0:
            data[11] = 1
            file_name = f"mouse_{data[0]}_session_{data[1]}_trial_{data[2]}.{data[3]}.v{data[4]}.{data[5]}.{data[6]}.{data[7]}.{data[8]}.{data[9]}.{data[10]}.{data[11]}"
            sql1 = "UPDATE Analysis SET registration_main=?,registration_v=? WHERE component_evaluation_main=? "
            val1 = [file_name, data[11], file[i]]
            cursor.execute(sql1, val1)

        else:
            data[11] += 1
            file_name = f"mouse_{data[0]}_session_{data[1]}_trial_{data[2]}.{data[3]}.v{data[4]}.{data[5]}.{data[6]}.{data[7]}.{data[8]}.{data[9]}.{data[10]}.{data[11]}"
            sql2 = "UPDATE Analysis SET registration_main=?,registration_v=?  WHERE motion_correction_main=? "
            val2 = [file_name, data[11], file[i]]
            cursor.execute(sql2, val2)
            database.commit()
        database.commit()

        output_file_path = data_dir + f'{file_name}.pkl'

    # Take alignment data for the timeline of alignment

        sql = "SELECT alignment_timeline FROM Analysis WHERE mouse = ? AND session=? AND component_evaluation_v=?"
        val = [mouse, session, component_evaluation_v]
        cursor.execute(sql, val)
        result = cursor.fetchall()
        alignment_timeline = []
        inter = []
        for x in result:
            inter += x
        for y in inter:
            alignment_timeline.append(y)

        # Multiple list created to append the relevant information for the registration and creation of a unique time trace
        # matrix (cnm.estimates.A  and cnm.estimates.C ) both taken after component evaluation

        A_list = []  # List for contour matrix on multiple trials
        FOV_size = []  # List for the cn filter dim (to verify it is always the same dims)
        A_number_components = []  # List with the total number of components extracted for each trial
        C_dims = []  # Dimension of C, to keep track of timeline
        C_list = []  # List with traces for each trial
        evaluated_trials = []
        typical_size = []

        if type(file[i]) == str:
            component_evaluation_hdf5_file_path = os.environ['DATA_DIR'] + '/interim/component_evaluation/trial_wise/' + \
                                                  file[i]
            sql = "SELECT source_extraction_corr FROM Analysis WHERE component_evaluation_main=?"
            val = [file[i]]
            cursor.execute(sql, val)
            result = cursor.fetchall()
            source_extraction_corr = []
            for x in result:
                source_extraction_corr += x
            corr_path = source_extraction_corr[i]
            cnm = load_CNMF(component_evaluation_hdf5_file_path)
            cn_filter = np.load(corr_path)
            FOV_size.append(cn_filter.shape)

            A_number_components.append(cnm.estimates.idx_components.shape[0])
            A_list.append(cnm.estimates.A[:, cnm.estimates.idx_components])
            C_dims.append(cnm.estimates.C.shape)
            size = cnm.estimates.A[:, cnm.estimates.idx_components].sum(axis=0)
            sql = "SELECT trial FROM Analysis WHERE component_evaluation_main=?"
            val = [file[i]]
            cursor.execute(sql, val)
            result = cursor.fetchall()
            trial=[]
            for x in result:
                trial += x
            sql = "SELECT is_rest FROM Analysis WHERE component_evaluation_main=?"
            val = [file[i]]
            cursor.execute(sql, val)
            result = cursor.fetchall()
            is_rest = []
            for x in result:
                is_rest += x
            evaluated_trials.append((trial[i] - 1) * 2 + is_rest[i])  # Number that goes from 0 to 42
            for j in range(len(cnm.estimates.idx_components)):
                typical_size.append(size[0, j])
            sql = "SELECT normalization FROM Analysis WHERE component_evaluation_main=?"
            val = [file[i]]
            cursor.execute(sql, val)
            result = cursor.fetchall()
            normalization_v = []
            for x in result:
                normalization_v += x
            if normalization_v[i] == 'True':
                if cnm.estimates.bl is None:
                    raw_normed, cnm_normed, res_normed, s_normed, noise_levels = normalization.normalize_traces(
                        cnm.estimates.C,
                        cnm.estimates.YrA,
                        cnm.estimates.S,
                        1,
                        offset_method="denoised_floor")
                else:
                    raw_normed, cnm_normed, res_normed, s_normed, noise_levels = normalization.normalize_traces(
                        cnm.estimates.C - cnm.estimates.bl[:, np.newaxis],
                        cnm.estimates.YrA,
                        cnm.estimates.S,
                        1,
                        offset_method="denoised_floor")
                C_list.append(cnm_normed[cnm.estimates.idx_components, :])
            else:
                if cnm.estimates.bl is None:
                    C_list.append(cnm.estimates.C[cnm.estimates.idx_components, :])
                else:
                    C_list.append(cnm.estimates.C[cnm.estimates.idx_components, :] - cnm.estimates.bl[
                        cnm.estimates.idx_components, np.newaxis])

        # Open the timeline and create the new traces matrix C_matrix
        with open(alignment_timeline[i], 'rb') as f:
            timeline = pickle.load(f)
        total_time1 = 0
        for i in range(len(C_list) - 1):
            total_time1 = total_time1 + C_list[i].shape[1]
        total_time2 = timeline[len(timeline) - 1][1] + C_list[np.argmax(evaluated_trials)].shape[1]
        total_time = max(total_time1, total_time2)
        timeline.append(['End', total_time])

        # Add a size restriction on the neurons that will further be processed. This restriction boundary
        # decision is based in the histogram of typical neuronal sizes
        new_A_list = []
        new_C_list = []
        A_components = []
        C_dims_new = []
        new_evaluated_trials = []
        for i in range(len(A_list)):
            accepted_size = []
            size = A_list[i].sum(axis=0)
            for j in range(size.shape[1]):
                if 10 < size[0, j] < 25:
                    accepted_size.append(j)
            if len(accepted_size) > 1:
                new_A_list.append(A_list[i][:, accepted_size])
                new_C_list.append(C_list[i][accepted_size, :])
                A_components.append(A_number_components[i])
                C_dims_new.append(new_C_list[-1].shape)
                new_evaluated_trials.append(evaluated_trials[i])
        A_list = new_A_list
        C_list = new_C_list

        # Run CaImAn registration routine that use the Hungarian matching algorithm in the contours list
        sql = "SELECT cost_threshold FROM Analysis WHERE mouse = ? AND session=? AND component_evaluation_v=?"
        val = [mouse, session, component_evaluation_v]
        cursor.execute(sql, val)
        result = cursor.fetchall()
        for x in result:
            cost_threshold = x
        sql = "SELECT max_dist FROM Analysis WHERE mouse = ? AND session=? AND component_evaluation_v=?"
        val = [mouse, session, component_evaluation_v]
        cursor.execute(sql, val)
        result = cursor.fetchall()
        for x in result:
            max_dist = x
        spatial_union, assignments, match = register_multisession(A=A_list, dims=FOV_size[0],
                                                                  thresh_cost=cost_threshold,
                                                                  max_dist=max_dist)

        C_matrix = np.zeros((spatial_union.shape[1], total_time))

        new_assignments = np.zeros((spatial_union.shape[1], len(timeline)))
        for i in range(spatial_union.shape[1]):
            for j in range(assignments.shape[1]):
                trial = new_evaluated_trials[j]
                if not math.isnan(assignments[i, j]):
                    new_assignments[i][trial] = assignments[i, j] + 1
        for i in range(spatial_union.shape[1]):
            for j in range(assignments.shape[1]):
                trial = new_evaluated_trials[j]
                if not math.isnan(assignments[i, j]):
                    C_matrix[i][timeline[trial][1]:timeline[trial][1] + (C_list[j])[int(assignments[i, j]), :].shape[0]] = (
                                                                                                                           C_list[
                                                                                                                               j])[
                                                                                                                           int(
                                                                                                                               assignments[
                                                                                                                                   i, j]),
                                                                                                                           :]

        cnm_registration = estimates(A=spatial_union, C=C_matrix)
        with open(output_file_path, 'wb') as output_file:
            pickle.dump(cnm_registration, output_file, pickle.HIGHEST_PROTOCOL)

    return
