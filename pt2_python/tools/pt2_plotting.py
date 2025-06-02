import os
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


plot_recipes = [
    {
        "name": "used_ls_pt",
        "title": "Used LS pT",
        "xrange": [0,250],
        "objects": ["ls_used_real_pt","ls_used_isfake_pt"],
        "legend": ["Real","Fake"],
        "xaxis_log":True,
        "nbins":180,
        "xlabel":"pT [GeV]",
        "ylabel":"Counts",
        "has_ratio":True,
        "ratio_log":True,
    },
    {
        "name": "used_ls_phi",
        "title": "Used LS phi",
        "xrange": [-3.2,3.2],
        "objects": ["ls_used_real_phi","ls_used_isfake_phi"],
        "legend": ["Real","Fake"],
        "xaxis_log":False,
        "nbins":180,
        "xlabel":"phi",
        "ylabel":"Counts",
        "has_ratio":True,
        "ratio_log":True,
    },
    {
        "name": "used_ls_eta",
        "title": "Used LS eta",
        "xrange": [-5,5],
        "objects": ["ls_used_real_eta","ls_used_isfake_eta"],
        "legend": ["Real","Fake"],
        "xaxis_log":False,
        "nbins":180,
        "xlabel":"eta",
        "ylabel":"Counts",
        "has_ratio":True,
        "ratio_log":True,
    },
    {
        "name": "unused_ls_pt",
        "title": "Unused LS pT",
        "xrange": [0,250],
        "objects": ["ls_unused_real_pt","ls_unused_isfake_pt"],
        "legend": ["Real","Fake"],
        "xaxis_log":True,
        "nbins":180,
        "xlabel":"pT [GeV]",
        "ylabel":"Counts",
        "has_ratio":True,
        "ratio_log":True,
    },
    {
        "name": "unused_ls_phi",
        "title": "Unused LS phi",
        "xrange": [-3.2,3.2],
        "objects": ["ls_unused_real_phi","ls_unused_isfake_phi"],
        "legend": ["Real","Fake"],
        "xaxis_log":False,
        "nbins":180,
        "xlabel":"phi",
        "ylabel":"Counts",
        "has_ratio":True,
        "ratio_log":True,
    },
    {
        "name": "unused_ls_eta",
        "title": "Unused LS eta",
        "xrange": [-5,5],
        "objects": ["ls_unused_real_eta","ls_unused_isfake_eta"],
        "legend": ["Real","Fake"],
        "xaxis_log":False,
        "nbins":180,
        "xlabel":"eta",
        "ylabel":"Counts",
        "has_ratio":True,
        "ratio_log":True,
    },
    {
        "name": "used_pls_pt",
        "title": "Used pLS pT",
        "xrange": [0,250],
        "objects": ["pls_used_real_pt","pls_used_isfake_pt"],
        "legend": ["Real","Fake"],
        "xaxis_log":True,
        "nbins":180,
        "xlabel":"pT [GeV]",
        "ylabel":"Counts",
        "has_ratio":True,
        "ratio_log":True,
    },
    {
        "name": "used_pls_phi",
        "title": "Used pLS phi",
        "xrange": [-3.2,3.2],
        "objects": ["pls_used_real_phi","pls_used_isfake_phi"],
        "legend": ["Real","Fake"],
        "xaxis_log":False,
        "nbins":180,
        "xlabel":"phi",
        "ylabel":"Counts",
        "has_ratio":True,
        "ratio_log":True,
    },
    {
        "name": "used_pls_eta",
        "title": "Used pLS eta",
        "xrange": [-5,5],
        "objects": ["pls_used_real_eta","pls_used_isfake_eta"],
        "legend": ["Real","Fake"],
        "xaxis_log":False,
        "nbins":180,
        "xlabel":"eta",
        "ylabel":"Counts",
        "has_ratio":True,
        "ratio_log":True,
    },
    {
        "name": "unused_pls_pt",
        "title": "Unused pLS pT",
        "xrange": [0,250],
        "objects": ["pls_unused_real_pt","pls_unused_isfake_pt"],
        "legend": ["Real","Fake"],
        "xaxis_log":True,
        "nbins":180,
        "xlabel":"pT [GeV]",
        "ylabel":"Counts",
        "has_ratio":True,
        "ratio_log":True,
    },
    {
        "name": "unused_pls_phi",
        "title": "Unused pLS phi",
        "xrange": [-3.2,3.2],
        "objects": ["pls_unused_real_phi","pls_unused_isfake_phi"],
        "legend": ["Real","Fake"],
        "xaxis_log":False,
        "nbins":180,
        "xlabel":"phi",
        "ylabel":"Counts",
        "has_ratio":True,
        "ratio_log":True,
    },
    {
        "name": "unused_pls_eta",
        "title": "Unused pLS eta",
        "xrange": [-5,5],
        "objects": ["pls_unused_real_eta","pls_unused_isfake_eta"],
        "legend": ["Real","Fake"],
        "xaxis_log":False,
        "nbins":180,
        "xlabel":"eta",
        "ylabel":"Counts",
        "has_ratio":True,
        "ratio_log":True,
    },
    {
        "name": "all_pls_pt",
        "title": "All pLS pT",
        "xrange": [0,250],
        "objects": ["pls_all_real_pt","pls_all_isfake_pt"],
        "legend": ["Real","Fake"],
        "xaxis_log":True,
        "nbins":180,
        "xlabel":"pT [GeV]",
        "ylabel":"Counts",
        "has_ratio":True,
        "ratio_log":True,
    },
    {
        "name": "all_pls_phi",
        "title": "All pLS phi",
        "xrange": [-3.2,3.2],
        "objects": ["pls_all_real_phi","pls_all_isfake_phi"],
        "legend": ["Real","Fake"],
        "xaxis_log":False,
        "nbins":180,
        "xlabel":"phi",
        "ylabel":"Counts",
        "has_ratio":True,
        "ratio_log":True,
    },
    {
        "name": "all_pls_eta",
        "title": "All pLS eta",
        "xrange": [-5,5],
        "objects": ["pls_all_real_eta","pls_all_isfake_eta"],
        "legend": ["Real","Fake"],
        "xaxis_log":False,
        "nbins":180,
        "xlabel":"eta",
        "ylabel":"Counts",
        "has_ratio":True,
        "ratio_log":True,
    },
    {
        "name": "all_ls_pt",
        "title": "All LS pT",
        "xrange": [0,250],
        "objects": ["ls_all_real_pt","ls_all_isfake_pt"],
        "legend": ["Real","Fake"],
        "xaxis_log":True,
        "nbins":180,
        "xlabel":"pT [GeV]",
        "ylabel":"Counts",
        "has_ratio":True,
        "ratio_log":True,
    },
    {
        "name": "all_ls_phi",
        "title": "All LS phi",
        "xrange": [-3.2,3.2],
        "objects": ["ls_all_real_phi","ls_all_isfake_phi"],
        "legend": ["Real","Fake"],
        "xaxis_log":False,
        "nbins":180,
        "xlabel":"phi",
        "ylabel":"Counts",
        "has_ratio":True,
        "ratio_log":True,
    },
    {
        "name": "all_ls_eta",
        "title": "All LS eta",
        "xrange": [-5,5],
        "objects": ["ls_all_real_eta","ls_all_isfake_eta"],
        "legend": ["Real","Fake"],
        "xaxis_log":False,
        "nbins":180,
        "xlabel":"eta",
        "ylabel":"Counts",
        "has_ratio":True,
        "ratio_log":True,
    },
    {
        "name": "used_ls_simIdx_breakdown_pt",
        "title": "Used LS pT SimIdx Matching",
        "xrange": [0,250],
        "objects": [
            "ls_used_real_simIdxmacthes_pls_unused_real_pt",
            "ls_used_real_simIdxmacthes_pls_used_real_pt",
            "ls_used_real_simIdxmacthes_none_pt",
            "ls_used_isfake_pt"
        ],
        "legend": [
            "Real and Matches Unused Real pLS",
            "Real and Matches Used Real pLS",
            "Real and No pLS Matches",
            "Fake"
        ],
        "xaxis_log":True,
        "nbins":180,
        "xlabel":"pT [GeV]",
        "ylabel":"Counts",
        "has_ratio":False,
        "ratio_log":False,
    },
    {
        "name": "used_ls_simIdx_breakdown_phi",
        "title": "Used LS phi SimIdx Matching",
        "xrange": [-3.2,3.2],
        "objects": [
            "ls_used_real_simIdxmacthes_pls_unused_real_phi",
            "ls_used_real_simIdxmacthes_pls_used_real_phi",
            "ls_used_real_simIdxmacthes_none_phi",
            "ls_used_isfake_phi"
        ],
        "legend": [
            "Real and Matches Unused Real pLS",
            "Real and Matches Used Real pLS",
            "Real and No pLS Matches",
            "Fake"
        ],
        "xaxis_log":False,
        "nbins":180,
        "xlabel":"phi",
        "ylabel":"Counts",
        "has_ratio":False,
        "ratio_log":False,
    },
    {
        "name": "used_ls_simIdx_breakdown_eta",
        "title": "Used LS eta SimIdx Matching",
        "xrange": [-5,5],
        "objects": [
            "ls_used_real_simIdxmacthes_pls_unused_real_eta",
            "ls_used_real_simIdxmacthes_pls_used_real_eta",
            "ls_used_real_simIdxmacthes_none_eta",
            "ls_used_isfake_eta"
        ],
        "legend": [
            "Real and Matches Unused Real pLS",
            "Real and Matches Used Real pLS",
            "Real and No pLS Matches",
            "Fake"
        ],
        "xaxis_log":False,
        "nbins":180,
        "xlabel":"eta",
        "ylabel":"Counts",
        "has_ratio":False,
        "ratio_log":False,
    },
    {
        "name": "unused_ls_simIdx_breakdown_pt",
        "title": "Unused LS pT SimIdx Matching",
        "xrange": [0,250],
        "objects": [
            "ls_unused_real_simIdxmacthes_pls_unused_real_pt",
            "ls_unused_real_simIdxmacthes_pls_used_real_pt",
            "ls_unused_real_simIdxmacthes_none_pt",
            "ls_unused_isfake_pt"
        ],
        "legend": [
            "Real and Matches Unused Real pLS",
            "Real and Matches Used Real pLS",
            "Real and No pLS Matches",
            "Fake"
        ],
        "xaxis_log":True,
        "nbins":180,
        "xlabel":"pT [GeV]",
        "ylabel":"Counts",
        "has_ratio":False,
        "ratio_log":False,
    },
    {
        "name": "unused_ls_simIdx_breakdown_phi",
        "title": "Unused LS phi SimIdx Matching",
        "xrange": [-3.2,3.2],
        "objects": [
            "ls_unused_real_simIdxmacthes_pls_unused_real_phi",
            "ls_unused_real_simIdxmacthes_pls_used_real_phi",
            "ls_unused_real_simIdxmacthes_none_phi",
            "ls_unused_isfake_phi"
        ],
        "legend": [
            "Real and Matches Unused Real pLS",
            "Real and Matches Used Real pLS",
            "Real and No pLS Matches",
            "Fake"
        ],
        "xaxis_log":False,
        "nbins":180,
        "xlabel":"phi",
        "ylabel":"Counts",
        "has_ratio":False,
        "ratio_log":False,
    },
    {
        "name": "unused_ls_simIdx_breakdown_eta",
        "title": "Unused LS eta SimIdx Matching",
        "xrange": [-5,5],
        "objects": [
            "ls_unused_real_simIdxmacthes_pls_unused_real_eta",
            "ls_unused_real_simIdxmacthes_pls_used_real_eta",
            "ls_unused_real_simIdxmacthes_none_eta",
            "ls_unused_isfake_eta"
        ],
        "legend": [
            "Real and Matches Unused Real pLS",
            "Real and Matches Used Real pLS",
            "Real and No pLS Matches",
            "Fake"
        ],
        "xaxis_log":False,
        "nbins":180,
        "xlabel":"eta",
        "ylabel":"Counts",
        "has_ratio":False,
        "ratio_log":False,
    },
    {
        "name": "used_pls_simIdx_breakdown_pt",
        "title": "Used pLS pT SimIdx Matching",
        "xrange": [0,250],
        "objects": [
            "pls_used_real_simIdxmacthes_ls_unused_real_pt",
            "pls_used_real_simIdxmacthes_ls_used_real_pt",
            "pls_used_real_simIdxmacthes_none_pt",
            "pls_used_isfake_pt"
        ],
        "legend": [
            "Real and Matches Unused Real LS",
            "Real and Matches Used Real LS",
            "Real and No LS Matches",
            "Fake"
        ],
        "xaxis_log":True,
        "nbins":180,
        "xlabel":"pT [GeV]",
        "ylabel":"Counts",
        "has_ratio":False,
        "ratio_log":False,
    },
    {
        "name": "used_pls_simIdx_breakdown_phi",
        "title": "Used pLS phi SimIdx Matching",
        "xrange": [-3.2,3.2],
        "objects": [
            "pls_used_real_simIdxmacthes_ls_unused_real_phi",
            "pls_used_real_simIdxmacthes_ls_used_real_phi",
            "pls_used_real_simIdxmacthes_none_phi",
            "pls_used_isfake_phi"
        ],
        "legend": [
            "Real and Matches Unused Real LS",
            "Real and Matches Used Real LS",
            "Real and No LS Matches",
            "Fake"
        ],
        "xaxis_log":False,
        "nbins":180,
        "xlabel":"phi",
        "ylabel":"Counts",
        "has_ratio":False,
        "ratio_log":False,
    },
    {
        "name": "used_pls_simIdx_breakdown_eta",
        "title": "Used pLS eta SimIdx Matching",
        "xrange": [-5,5],
        "objects": [
            "pls_used_real_simIdxmacthes_ls_unused_real_eta",
            "pls_used_real_simIdxmacthes_ls_used_real_eta",
            "pls_used_real_simIdxmacthes_none_eta",
            "pls_used_isfake_eta"
        ],
        "legend": [
            "Real and Matches Unused Real LS",
            "Real and Matches Used Real LS",
            "Real and No LS Matches",
            "Fake"
        ],
        "xaxis_log":False,
        "nbins":180,
        "xlabel":"eta",
        "ylabel":"Counts",
        "has_ratio":False,
        "ratio_log":False,
    },
    {
        "name": "unused_pls_simIdx_breakdown_pt",
        "title": "Unused pLS pT SimIdx Matching",
        "xrange": [0,250],
        "objects": [
            "pls_unused_real_simIdxmacthes_ls_unused_real_pt",
            "pls_unused_real_simIdxmacthes_ls_used_real_pt",
            "pls_unused_real_simIdxmacthes_none_pt",
            "pls_unused_isfake_pt"
        ],
        "legend": [
            "Real and Matches Unused Real LS",
            "Real and Matches Used Real LS",
            "Real and No LS Matches",
            "Fake"
        ],
        "xaxis_log":True,
        "nbins":180,
        "xlabel":"pT [GeV]",
        "ylabel":"Counts",
        "has_ratio":False,
        "ratio_log":False,
    },
    {
        "name": "unused_pls_simIdx_breakdown_phi",
        "title": "Unused pLS phi SimIdx Matching",
        "xrange": [-3.2,3.2],
        "objects": [
            "pls_unused_real_simIdxmacthes_ls_unused_real_phi",
            "pls_unused_real_simIdxmacthes_ls_used_real_phi",
            "pls_unused_real_simIdxmacthes_none_phi",
            "pls_unused_isfake_phi"
        ],
        "legend": [
            "Real and Matches Unused Real LS",
            "Real and Matches Used Real LS",
            "Real and No LS Matches",
            "Fake"
        ],
        "xaxis_log":False,
        "nbins":180,
        "xlabel":"phi",
        "ylabel":"Counts",
        "has_ratio":False,
        "ratio_log":False,
    },
    {
        "name": "unused_pls_simIdx_breakdown_eta",
        "title": "Unused pLS eta SimIdx Matching",
        "xrange": [-5,5],
        "objects": [
            "pls_unused_real_simIdxmacthes_ls_unused_real_eta",
            "pls_unused_real_simIdxmacthes_ls_used_real_eta",
            "pls_unused_real_simIdxmacthes_none_eta",
            "pls_unused_isfake_eta"
        ],
        "legend": [
            "Real and Matches Unused Real LS",
            "Real and Matches Used Real LS",
            "Real and No LS Matches",
            "Fake"
        ],
        "xaxis_log":False,
        "nbins":180,
        "xlabel":"eta",
        "ylabel":"Counts",
        "has_ratio":False,
        "ratio_log":False,
    },
    {
        "name": "used_ls_simIdx_breakdown_pt_log",
        "title": "Used LS pT SimIdx Matching",
        "xrange": [0,250],
        "objects": [
            "ls_used_real_simIdxmacthes_pls_unused_real_pt",
            "ls_used_real_simIdxmacthes_pls_used_real_pt",
            "ls_used_real_simIdxmacthes_none_pt",
            "ls_used_isfake_pt"
        ],
        "legend": [
            "Real and Matches Unused Real pLS",
            "Real and Matches Used Real pLS",
            "Real and No pLS Matches",
            "Fake"
        ],
        "xaxis_log":True,
        "yaxis_log":True,
        "nbins":180,
        "xlabel":"pT [GeV]",
        "ylabel":"Counts",
        "has_ratio":False,
        "ratio_log":False,
    },
    {
        "name": "used_ls_simIdx_breakdown_phi_log",
        "title": "Used LS phi SimIdx Matching",
        "xrange": [-3.2,3.2],
        "objects": [
            "ls_used_real_simIdxmacthes_pls_unused_real_phi",
            "ls_used_real_simIdxmacthes_pls_used_real_phi",
            "ls_used_real_simIdxmacthes_none_phi",
            "ls_used_isfake_phi"
        ],
        "legend": [
            "Real and Matches Unused Real pLS",
            "Real and Matches Used Real pLS",
            "Real and No pLS Matches",
            "Fake"
        ],
        "xaxis_log":False,
        "yaxis_log":True,
        "nbins":180,
        "xlabel":"phi",
        "ylabel":"Counts",
        "has_ratio":False,
        "ratio_log":False,
    },
    {
        "name": "used_ls_simIdx_breakdown_eta_log",
        "title": "Used LS eta SimIdx Matching",
        "xrange": [-5,5],
        "objects": [
            "ls_used_real_simIdxmacthes_pls_unused_real_eta",
            "ls_used_real_simIdxmacthes_pls_used_real_eta",
            "ls_used_real_simIdxmacthes_none_eta",
            "ls_used_isfake_eta"
        ],
        "legend": [
            "Real and Matches Unused Real pLS",
            "Real and Matches Used Real pLS",
            "Real and No pLS Matches",
            "Fake"
        ],
        "xaxis_log":False,
        "yaxis_log":True,
        "nbins":180,
        "xlabel":"eta",
        "ylabel":"Counts",
        "has_ratio":False,
        "ratio_log":False,
    },

    {
        "name": "unused_ls_simIdx_breakdown_pt_log",
        "title": "Unused LS pT SimIdx Matching",
        "xrange": [0,250],
        "objects": [
            "ls_unused_real_simIdxmacthes_pls_unused_real_pt",
            "ls_unused_real_simIdxmacthes_pls_used_real_pt",
            "ls_unused_real_simIdxmacthes_none_pt",
            "ls_unused_isfake_pt"
        ],
        "legend": [
            "Real and Matches Unused Real pLS",
            "Real and Matches Used Real pLS",
            "Real and No pLS Matches",
            "Fake"
        ],
        "xaxis_log":True,
        "yaxis_log":True,
        "nbins":180,
        "xlabel":"pT [GeV]",
        "ylabel":"Counts",
        "has_ratio":False,
        "ratio_log":False,
    },
    {
        "name": "unused_ls_simIdx_breakdown_phi_log",
        "title": "Unused LS phi SimIdx Matching",
        "xrange": [-3.2,3.2],
        "objects": [
            "ls_unused_real_simIdxmacthes_pls_unused_real_phi",
            "ls_unused_real_simIdxmacthes_pls_used_real_phi",
            "ls_unused_real_simIdxmacthes_none_phi",
            "ls_unused_isfake_phi"
        ],
        "legend": [
            "Real and Matches Unused Real pLS",
            "Real and Matches Used Real pLS",
            "Real and No pLS Matches",
            "Fake"
        ],
        "xaxis_log":False,
        "yaxis_log":True,
        "nbins":180,
        "xlabel":"phi",
        "ylabel":"Counts",
        "has_ratio":False,
        "ratio_log":False,
    },
    {
        "name": "unused_ls_simIdx_breakdown_eta_log",
        "title": "Unused LS eta SimIdx Matching",
        "xrange": [-5,5],
        "objects": [
            "ls_unused_real_simIdxmacthes_pls_unused_real_eta",
            "ls_unused_real_simIdxmacthes_pls_used_real_eta",
            "ls_unused_real_simIdxmacthes_none_eta",
            "ls_unused_isfake_eta"
        ],
        "legend": [
            "Real and Matches Unused Real pLS",
            "Real and Matches Used Real pLS",
            "Real and No pLS Matches",
            "Fake"
        ],
        "xaxis_log":False,
        "yaxis_log":True,
        "nbins":180,
        "xlabel":"eta",
        "ylabel":"Counts",
        "has_ratio":False,
        "ratio_log":False,
    },
    {
        "name": "used_pls_simIdx_breakdown_pt_log",
        "title": "Used pLS pT SimIdx Matching",
        "xrange": [0,250],
        "objects": [
            "pls_used_real_simIdxmacthes_ls_unused_real_pt",
            "pls_used_real_simIdxmacthes_ls_used_real_pt",
            "pls_used_real_simIdxmacthes_none_pt",
            "pls_used_isfake_pt"
        ],
        "legend": [
            "Real and Matches Unused Real LS",
            "Real and Matches Used Real LS",
            "Real and No LS Matches",
            "Fake"
        ],
        "xaxis_log":True,
        "yaxis_log":True,
        "nbins":180,
        "xlabel":"pT [GeV]",
        "ylabel":"Counts",
        "has_ratio":False,
        "ratio_log":False,
    },
    {
        "name": "used_pls_simIdx_breakdown_phi_log",
        "title": "Used pLS phi SimIdx Matching",
        "xrange": [-3.2,3.2],
        "objects": [
            "pls_used_real_simIdxmacthes_ls_unused_real_phi",
            "pls_used_real_simIdxmacthes_ls_used_real_phi",
            "pls_used_real_simIdxmacthes_none_phi",
            "pls_used_isfake_phi"
        ],
        "legend": [
            "Real and Matches Unused Real LS",
            "Real and Matches Used Real LS",
            "Real and No LS Matches",
            "Fake"
        ],
        "xaxis_log":False,
        "yaxis_log":True,
        "nbins":180,
        "xlabel":"phi",
        "ylabel":"Counts",
        "has_ratio":False,
        "ratio_log":False,
    },
    {
        "name": "used_pls_simIdx_breakdown_eta_log",
        "title": "Used pLS eta SimIdx Matching",
        "xrange": [-5,5],
        "objects": [
            "pls_used_real_simIdxmacthes_ls_unused_real_eta",
            "pls_used_real_simIdxmacthes_ls_used_real_eta",
            "pls_used_real_simIdxmacthes_none_eta",
            "pls_used_isfake_eta"
        ],
        "legend": [
            "Real and Matches Unused Real LS",
            "Real and Matches Used Real LS",
            "Real and No LS Matches",
            "Fake"
        ],
        "xaxis_log":False,
        "yaxis_log":True,
        "nbins":180,
        "xlabel":"eta",
        "ylabel":"Counts",
        "has_ratio":False,
        "ratio_log":False,
    },
    {
        "name": "unused_pls_simIdx_breakdown_pt_log",
        "title": "Unused pLS pT SimIdx Matching",
        "xrange": [0,250],
        "objects": [
            "pls_unused_real_simIdxmacthes_ls_unused_real_pt",
            "pls_unused_real_simIdxmacthes_ls_used_real_pt",
            "pls_unused_real_simIdxmacthes_none_pt",
            "pls_unused_isfake_pt"
        ],
        "legend": [
            "Real and Matches Unused Real LS",
            "Real and Matches Used Real LS",
            "Real and No LS Matches",
            "Fake"
        ],
        "xaxis_log":True,
        "yaxis_log":True,
        "nbins":180,
        "xlabel":"pT [GeV]",
        "ylabel":"Counts",
        "has_ratio":False,
        "ratio_log":False,
    },
    {
        "name": "unused_pls_simIdx_breakdown_phi_log",
        "title": "Unused pLS phi SimIdx Matching",
        "xrange": [-3.2,3.2],
        "objects": [
            "pls_unused_real_simIdxmacthes_ls_unused_real_phi",
            "pls_unused_real_simIdxmacthes_ls_used_real_phi",
            "pls_unused_real_simIdxmacthes_none_phi",
            "pls_unused_isfake_phi"
        ],
        "legend": [
            "Real and Matches Unused Real LS",
            "Real and Matches Used Real LS",
            "Real and No LS Matches",
            "Fake"
        ],
        "xaxis_log":False,
        "yaxis_log":True,
        "nbins":180,
        "xlabel":"phi",
        "ylabel":"Counts",
        "has_ratio":False,
        "ratio_log":False,
    },
    {
        "name": "unused_pls_simIdx_breakdown_eta_log",
        "title": "Unused pLS eta SimIdx Matching",
        "xrange": [-5,5],
        "objects": [
            "pls_unused_real_simIdxmacthes_ls_unused_real_eta",
            "pls_unused_real_simIdxmacthes_ls_used_real_eta",
            "pls_unused_real_simIdxmacthes_none_eta",
            "pls_unused_isfake_eta"
        ],
        "legend": [
            "Real and Matches Unused Real LS",
            "Real and Matches Used Real LS",
            "Real and No LS Matches",
            "Fake"
        ],
        "xaxis_log":False,
        "yaxis_log":False,
        "nbins":180,
        "xlabel":"eta",
        "ylabel":"Counts",
        "has_ratio":False,
        "ratio_log":False,
    },

]


def make_plots(data_dict, output_dir):
 
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
 
    # Define custom colors (colorblind friendly)
    colors = ['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
 
    for recipe in plot_recipes:
        name = recipe["name"]
        title = recipe["title"]
        objects = recipe["objects"]
        legend = recipe["legend"]
        nbins = recipe["nbins"]
        rangex = recipe["xrange"]
        xlabel = recipe["xlabel"]
        ylabel = recipe["ylabel"]
        xaxis_log = recipe["xaxis_log"]
        has_ratio = recipe["has_ratio"]
        ratio_log = recipe["ratio_log"]
        yaxis_log = recipe.get("yaxis_log", False)  # <-- NEW LINE
 
        if not has_ratio and ratio_log:
            raise Exception(f"Plot '{name}' has `ratio_log`=True but `has_ratio`=False.")
 
        # Pull and flatten the data
        hist_data = [ak.to_numpy(ak.flatten(data_dict[obj])) for obj in objects]
 
        # Bin edges for consistent binning
        if xaxis_log:
            bin_edges = np.logspace(np.log10(max(rangex[0],0.1)), np.log10(rangex[1]), nbins + 1)
        else:
            bin_edges = np.linspace(rangex[0], rangex[1], nbins + 1)
 
        # Compute histograms
        hist_values = [np.histogram(data, bins=bin_edges)[0] for data in hist_data]
 
        # Start figure
        if has_ratio:
            fig = plt.figure(figsize=(8, 6))
            gs = GridSpec(2, 1, height_ratios=[3, 1], hspace=0.05)
            ax_main = fig.add_subplot(gs[0])
            ax_ratio = fig.add_subplot(gs[1], sharex=ax_main)
        else:
            fig, ax_main = plt.subplots(figsize=(8, 5))
 
        # Plot stacked histogram
        ax_main.hist(
            hist_data,
            bins=bin_edges,
            stacked=True,
            label=legend,
            color=colors[:len(hist_data)],
            histtype='stepfilled',
            edgecolor='black',
        )
 
        ax_main.set_ylabel(ylabel)
        ax_main.set_title(title)
        ax_main.legend()
 
        if xaxis_log:
            ax_main.set_xscale('log')
            ax_main.axvline(x=0.8, color='red', linestyle='--', linewidth=1.5, label='x = 0.8')
 
        if yaxis_log:  # <-- NEW BLOCK
            ax_main.set_yscale('log')
 
        if has_ratio:
            # Compute ratio: first / second
            ratio = np.zeros_like(hist_values[0], dtype=float)
            with np.errstate(divide='ignore', invalid='ignore'):
                ratio = np.true_divide(hist_values[0], hist_values[1])
                ratio[hist_values[1] == 0] = np.nan
 
            bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
            ax_ratio.plot(bin_centers, ratio, drawstyle='steps-mid', color='black')
 
            ax_ratio.set_xlabel(xlabel)
            ax_ratio.set_ylabel('Ratio')
 
            if ratio_log:
                ax_ratio.set_yscale('log')
                ax_ratio.set_ylim(1e-3, 1.2)
                ax_ratio.set_yticks([1, 0.1, 0.01, 0.001])
                ax_ratio.get_yaxis().set_major_formatter(plt.ScalarFormatter())
            else:
                ax_ratio.set_ylim(0, 2)
 
        else:
            ax_main.set_xlabel(xlabel)
 
        # Save the plot
        output_path = os.path.join(output_dir, f"{name}.pdf")
        if has_ratio:
            fig.subplots_adjust(hspace=0.05)  # Already done by GridSpec, but reassert just in case
        else:
            fig.tight_layout()
 
        plt.tight_layout()
        plt.savefig(output_path)
 
        plt.close()
        print(f"Plot {name}.pdf has been created!")





#def make_plots(data_dict, output_dir):
#
#    # Ensure output directory exists
#    os.makedirs(output_dir, exist_ok=True)
#
#    # Define custom colors (colorblind friendly)
#    colors = ['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
#
#    for recipe in plot_recipes:
#        name = recipe["name"]
#        title = recipe["title"]
#        objects = recipe["objects"]
#        legend = recipe["legend"]
#        nbins = recipe["nbins"]
#        rangex = recipe["xrange"]
#        xlabel = recipe["xlabel"]
#        ylabel = recipe["ylabel"]
#        xaxis_log = recipe["xaxis_log"]
#        has_ratio = recipe["has_ratio"]
#        ratio_log = recipe["ratio_log"]
#
#        if not has_ratio and ratio_log:
#            raise Exception(f"Plot '{name}' has `ratio_log`=True but `has_ratio`=False.")
#
#        # Pull and flatten the data
#        hist_data = [ak.to_numpy(ak.flatten(data_dict[obj])) for obj in objects]
#
#        # Bin edges for consistent binning
#        if xaxis_log:
#            bin_edges = np.logspace(np.log10(max(rangex[0],0.1)), np.log10(rangex[1]), nbins + 1)
#        else:
#            bin_edges = np.linspace(rangex[0], rangex[1], nbins + 1)
#
#        # Compute histograms
#        hist_values = [np.histogram(data, bins=bin_edges)[0] for data in hist_data]
#
#        # Start figure
#        if has_ratio:
#            fig = plt.figure(figsize=(8, 6))
#            gs = GridSpec(2, 1, height_ratios=[3, 1], hspace=0.05)
#            ax_main = fig.add_subplot(gs[0])
#            ax_ratio = fig.add_subplot(gs[1], sharex=ax_main)
#        else:
#            fig, ax_main = plt.subplots(figsize=(8, 5))
#
#        # Plot stacked histogram
#        ax_main.hist(
#            hist_data,
#            bins=bin_edges,
#            stacked=True,
#            label=legend,
#            color=colors[:len(hist_data)],
#            histtype='stepfilled',
#            edgecolor='black',
#        )
#
#        ax_main.set_ylabel(ylabel)
#        ax_main.set_title(title)
#        ax_main.legend()
#
#        if xaxis_log:
#            ax_main.set_xscale('log')
#            ax_main.axvline(x=0.8, color='red', linestyle='--', linewidth=1.5, label='x = 0.8')
#
#        if has_ratio:
#            # Compute ratio: first / second
#            ratio = np.zeros_like(hist_values[0], dtype=float)
#            with np.errstate(divide='ignore', invalid='ignore'):
#                ratio = np.true_divide(hist_values[0], hist_values[1])
#                ratio[hist_values[1] == 0] = np.nan
#
#            bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
#            ax_ratio.plot(bin_centers, ratio, drawstyle='steps-mid', color='black')
#
#            ax_ratio.set_xlabel(xlabel)
#            ax_ratio.set_ylabel('Ratio')
#
#            if ratio_log:
#                ax_ratio.set_yscale('log')
#                ax_ratio.set_ylim(1e-3, 1.2)
#                ax_ratio.set_yticks([1, 0.1, 0.01, 0.001])
#                ax_ratio.get_yaxis().set_major_formatter(plt.ScalarFormatter())
#            else:
#                ax_ratio.set_ylim(0, 2)
#
#        else:
#            ax_main.set_xlabel(xlabel)
#
#        # Save the plot
#        output_path = os.path.join(output_dir, f"{name}.pdf")
#        if has_ratio:
#            fig.subplots_adjust(hspace=0.05)  # Already done by GridSpec, but reassert just in case
#        else:
#            fig.tight_layout()
#
#        plt.tight_layout()
#        plt.savefig(output_path)
#
#        plt.close()
#        print(f"Plot {name}.pdf has been created!")
