import os

from voyager_sdk.operator.operator import Operator
from voyager_sdk.file_repository import FileRepository
from voyager_sdk.protocols.processors.file_processor import FileProcessor


class TempoSomatic(Operator):
    """Base operator class for the SDK"""

    def __init__(
            self,
            request_id=None,
            run_ids=[],
            pipeline=None,
            pairing=None,
            output_directory_prefix=None,
            file_group=None,
            job_group_id=None,
            job_group_notifier_id=None,
            **kwargs):
        """
        request_id: metadata key:igoRequestId
        run_ids: ids of previous runs
        pipeline: {
            "pipeline_format": "",
            "pipeline_link": "",
            "pipeline_version: "",
            "pipeline_entrypoint": ""
        },
        file_group: file_group_id
        """
        self.tempo_file_group = os.environ.get("TEMPO_BAM")
        self.pairing_file_location = os.environ.get("TEMPO_PAIRING_FILE")
        super().__init__(
            request_id=request_id,
            run_ids=run_ids or [],
            pipeline=pipeline,
            pairing=pairing,
            output_directory_prefix=output_directory_prefix,
            file_group=file_group,
            job_group_id=job_group_id,
            job_group_notifier_id=job_group_notifier_id,
            **kwargs)
        
    def get_jobs(self):
        """
        :return:  A list of dicts where each dict contains
            :app: dict in pipeline format
            :inputs: dict of inputs for the pipeline run
            :name: run name (Example: TumorSampleName Run 1)
        """
        jobs = []
        for run in self.get_runs():
            cmo_sample_name = run.tags["cmoSampleName"]
            target = run.tags["target"]
            bam = next(FileRepository.filter(metadata=f"cmoSampleName:{cmo_sample_name}", file_type="bam", file_group=self.tempo_file_group), None)
            bai = next(FileRepository.filter(metadata=f"cmoSampleName:{cmo_sample_name}", file_type="bai", file_group=self.tempo_file_group), None)
            if bam and bai:
                tumorOrNormal = bam.metadata["tumorOrNormal"]
                if tumorOrNormal == "Tumor":
                    # Find paired normal
                    normal = self.find_normal(cmo_sample_name)
                    if normal:
                        normal_bam = next(FileRepository.filter(metadata=f"cmoSampleName:{normal}", file_type="bam",
                                                                file_group=self.tempo_file_group), None)
                        normal_bai = next(FileRepository.filter(metadata=f"cmoSampleName:{normal}", file_type="bai",
                                                                file_group=self.tempo_file_group), None)
                        job = self.create_job(target, cmo_sample_name, bam, bai, normal, normal_bam, normal_bai)
                        jobs.append(job)
                else:
                    # Find paired tumors
                    tumors = self.find_tumors(cmo_sample_name)
                    for tumor in tumors:
                        tumor_bam = FileRepository.filter(metadata=f"cmoSampleName:{tumor}", file_type="bam",
                                                          file_group=self.tempo_file_group)
                        tumor_bai = FileRepository.filter(metadata=f"cmoSampleName:{tumor}", file_type="bai",
                                                          file_group=self.tempo_file_group)
                        job = self.create_job(target, cmo_sample_name, tumor_bam, tumor_bai, cmo_sample_name, bam, bai)
                        jobs.append(job)
        return jobs


    def _get_pairs_from_file(self):
        """
        Parse the pairing file (header: NORMAL_ID TUMOR_ID, whitespace separated).
        :return: list of unique (normal, tumor) tuples, in file order
        """
        pairs = []
        seen = set()
        with open(self.pairing_file_location) as f:
            header = f.readline().split()
            if [h.upper() for h in header] != ["NORMAL_ID", "TUMOR_ID"]:
                raise ValueError(f"Unexpected pairing file header: {header}")
            for line in f:
                cols = line.split()
                if not cols:
                    continue
                if len(cols) != 2:
                    raise ValueError(f"Malformed pairing line: {line!r}")
                pair = (cols[0], cols[1])
                if pair not in seen:
                    seen.add(pair)
                    pairs.append(pair)
        return pairs

    def find_normal(self, tumor):
        """
        :return: the normal sample paired with the tumor, or None if not paired
        """
        for normal, t in self._get_pairs_from_file():
            if t == tumor:
                return normal
        return None

    def find_tumors(self, normal):
        """
        :return: list of tumor samples paired with the normal (empty if none)
        """
        return [t for n, t in self._get_pairs_from_file() if n == normal]

    def create_job(self, target, tumor, tumor_bam, tumor_bai, normal, normal_bam, normal_bai):
        inputs = {
            "somatic": True,
            "aggregate": False,
            "assayType": "exome",
            "genome": None,
            "workflows": "snv,sv,mutsig,lohhla,facets,msisensor,germSNV,germSV",
            "anonymizeFQ": True,
            "bamMapping": [
                {
                    "sample": tumor,
                    "target": target,
                    "bam": {
                        "class": "File",
                        "location": f"iris://{tumor_bam.path}"
                    },
                    "bai": {
                        "class": "File",
                        "location": f"iris://{tumor_bai.path}"
                    }
                },
                {
                    "sample": normal,
                    "target": target,
                    "bam": {
                        "class": "File",
                        "location": f"iris://{normal_bam.path}"
                    },
                    "bai": {
                        "class": "File",
                        "location": f"iris://{normal_bai.path}"
                    }
                }
            ],
            "pairing": [
                {
                    "tumor_id": tumor,
                    "normal_id": normal
                }
            ]
        }
        return inputs
