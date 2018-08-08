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
        self.barcode_id = str(run_info['barcodeId'])
        self.base_caller_args = str(run_info['basecallerArgs'])
        self.chip_type = str(run_info['chiptype'])
        self.experiment_name = str(run_info['expName'])
        self.bead_find_args = str(run_info['beadfindArgs'])
        self.do_base_recalibration = str(run_info['doBaseRecal'])
        self.site_name = str(run_info['site_name'])
        self.user_name = str(run_info['username'])
        self.tmap_version = str(run_info['tmap_version'])
        self.chef_chip_expiration_1 = str(run_info['exp_json']['chefChipExpiration1'])
        self.chef_chip_expiration_2 = str(run_info['exp_json']['chefChipExpiration2'])
        self.chef_chip_type_1 = str(run_info['exp_json']['chefChipType1'])
        self.chef_chip_type_2 = str(run_info['exp_json']['chefChipType2'])
        self.chef_instrument_name = str(run_info['exp_json']['chefInstrumentName'])
        self.chef_kit_type = str(run_info['exp_json']['chefKitType'])
        self.chef_last_update = datetime.strptime(str(run_info['exp_json']['chefLastUpdate']),
                                                  "%Y-%m-%dT%H:%M:%SZ").isoformat()
        self.chef_log_path = str(run_info['exp_json']['chefLogPath'])
        self.chef_lot_number = str(run_info['exp_json']['chefLotNumber'])
        self.chef_manufacture_date = str(run_info['exp_json']['chefManufactureDate'])
        self.chef_message = str(run_info['exp_json']['chefMessage'])
        self.chef_package_version = str(run_info['exp_json']['chefPackageVer'])
        self.chef_reagent_id = str(run_info['exp_json']['chefReagentID'])
        self.chef_reagents_expiration = str(run_info['exp_json']['chefReagentsExpiration'])
        self.chef_reagents_lot = str(run_info['exp_json']['chefReagentsLot'])
        self.chef_reagents_part = str(run_info['exp_json']['chefReagentsPart'])
        self.chef_reagents_serial_number = str(run_info['exp_json']['chefReagentsSerialNum'])
        self.chef_script_version = str(run_info['exp_json']['chefScriptVersion'])
        self.chef_solutions_expiration = str(run_info['exp_json']['chefSolutionsExpiration'])
        self.chef_solutions_lot = str(run_info['exp_json']['chefSolutionsLot'])
        self.chef_solutions_part = str(run_info['exp_json']['chefSolutionsPart'])
        self.chef_solutions_serial_number = str(run_info['exp_json']['chefSolutionsSerialNum'])
        self.chef_status = str(run_info['exp_json']['chefStatus'])
        self.chef_tip_rack_barcode = str(run_info['exp_json']['chefTipRackBarcode'])
        self.chip_barcode = str(run_info['exp_json']['chipBarcode'])
        self.chip_type = str(run_info['exp_json']['chipType'])
        self.pcr_cycles = str(run_info['exp_json']['cycles'])
        self.date = datetime.strptime(str(run_info['exp_json']['date']), "%Y-%m-%dT%H:%M:%SZ").isoformat()
        self.display_name = str(run_info['exp_json']['displayName'])
        self.experiment_directory = str(run_info['exp_json']['expDir'])
        self.flows = str(run_info['exp_json']['flows'])
        self.server = str(run_info['net_location'])
        self.run_id = str(run_info['exp_json']['id'])
        self.pgm_name = str(run_info['exp_json']['pgmName'])

        self.platform = str(run_info['exp_json']['platform'])
        self.s5_version = str(run_info['exp_json']['log']['s5_release_version'])
        self.reagent_barcode = str(run_info['exp_json']['reagentBarcode'])
        self.run_mode = str(run_info['exp_json']['runMode'])
        self.sequence_kit_barcode = str(run_info['exp_json']['seqKitBarcode'])
        self.sequence_kit_name = str(run_info['exp_json']['sequencekitname'])
        self.templating_size = str(run_info['plan']['templatingSize'])
        self.library_read_length = str(run_info['plan']['libraryReadLength'])
        self.transferred_run = "runTransferFromSource" in str(run_info['plan']['metaData'])
        self.analysis_args = str(run_info['experimentAnalysisSettings']['analysisargs'])
        self.base_recalibration_mode = str(run_info['experimentAnalysisSettings']['base_recalibration_mode'])
        self.base_recalibration_args = str(run_info['experimentAnalysisSettings']['calibrateargs'])
        self.library_kit_name = str(run_info['experimentAnalysisSettings']['libraryKitName'])

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
