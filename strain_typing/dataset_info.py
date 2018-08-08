import json
import os


class BaseCallerStats(object):

    def __init__(self, run_directory):
        self.basecaller_info = json.load(open(os.path.join(run_directory, "basecaller_results",
                                                           "datasets_basecaller.json")))

        self.barcode_stats = {}
        for k, v in self.basecaller_info['read_groups'].items():
            self.barcode_stats[k.split(".")[-1]] = {"total_bases": int(v["total_bases"]),
                                                    "read_count": int(v["read_count"]),
                                                    "Q20_bases": int(v["Q20_bases"])}
