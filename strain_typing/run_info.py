import json
import os
import sys
from datetime import datetime
from pprint import pprint


class RunInfo(object):
    """
    The purpose of this class is to hold the run/plan level information so that checks can be performed
    It also helps me to see what is available for QA
    """
    def __init__(self, run_directory):
        # ion_params_00.json
        run_info = json.load(open(os.path.join(run_directory, "ion_params_00.json")))
        self.date = datetime.strptime(str(run_info['exp_json']['date']), "%Y-%m-%dT%H:%M:%SZ").isoformat()

        self._warnings = None
        if "__alarm" in run_info['exp_json']['log']:
            self._warnings = str(run_info['exp_json']['log']['__alarm'])

        self.chip_log = run_info['exp_json']['log']

    @property
    def warnings(self):
        return self._warnings

    def check_run_info(self):
        """
        runs a serious of value check to make sure the run plan has not been modified from validated settings
        :yields: passed (bool), expected_result (str), actual_result (str), test_name (str)
        """
        tests_to_perform = {
            # "chef_kit_type": "Ion 510 & Ion 520 & Ion 530 Kit-Chef",
            # "sequence_kit_name": "Ion S5 Sequencing Kit",
            # 'library_kit_name': "Ion Xpress Plus Fragment Library Kit",
            # "library_read_length": "400",
            "templating_size": "400",
            "barcode_id": "IonCode",
            "chip_type": "530",
            "flows": "850",
            "base_caller_args": "BaseCaller --barcode-filter 0.01 --barcode-filter-minreads 10 "
                                "--phasing-residual-filter=2.0 --num-unfiltered 1000 "
                                "--barcode-filter-postpone 1 --qual-filter true --qual-filter-slope 0.040 "
                                "--qual-filter-offset 1.0 --wells-normalization on",
            "bead_find_args": "justBeadFind --args-json /opt/ion/config/args_530_beadfind.json",
            "do_base_recalibration": "panel_recal",
            "base_recalibration_mode": "panel_recal",
            "base_recalibration_args": "Calibration --num-calibration-regions 1,1",
            's5_version': "5.8"
        }

        test_to_skip_when_transferred = ["library_read_length", "templating_size"]
        for attribute, expected_result in tests_to_perform.items():
            if attribute in test_to_skip_when_transferred and self.transferred_run:
                sys.stderr.write("WARN: SKIPPING '{0}' TEST BECAUSE RUN HAS BEEN TRANSFERRED\n".format(attribute))
                continue
            yield self.__check_attr(attribute, expected_result)

    def __check_attr(self, attribute, expected_result):
        """
        private method used by 'check_run_info'
        :param attribute: the name of the attribute to test
        :param expected_result: the expected value of the attribute
        :return: passed (bool), expected_result (str), actual_result (str), test_name (str)
        """
        return getattr(self, attribute) == expected_result, expected_result, getattr(self, attribute), attribute
