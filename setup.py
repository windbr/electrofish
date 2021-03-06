__author__ = 'jpresern'
from distutils.core import setup

setup(name='efish',
      version='1.00',
      author="Janez Presern",
      description="Python 3 scripts for analysis of RELACS obtained neurophys data",
      keywords="neuro",
      # packages=['efish'],
      py_module=['efish'],
      # py_modules=['coherence',
      #             'eFish_analysis',
      #             'FI',
      #             'html_dump',
      #             'markup',
      #             'mathEngine',
      #             'parallel_analysis',
      #             'parse_voltage_trace_from_txt'
      #             'read_info_dat',
      #             'report_coherence_dyn',
      #             'report_FI',
      #             'report_FI_curves_dyn',
      #             'report_FI_freq',
      #             'report_FI_trace_dyn',
      #             'report_noise_traces',
      #             'report_resistance',
      #             'report_whitenoise_trace_dyn',
      #             'resistance',
      #             'single_analysis',
      #             'spike_extraction',
      #             'spont_act',
      #             'transfer_curve',
      #             'utils',
      #             'VI'],
      requires=['matplotlib', 'numpy', 'scipy', 'pyrelacs', 'PeakUtils']
      )