from collections import UserDict
import json

import logging as log


class Sample:
    # key: column name; value: column index in a ped9 tab-separated file according to spec
    ATTRIBUTES = {
        "family_id": 0,
        "individual_id": 1,
        "paternal_id": 2,
        "maternal_id": 3,
        "sex": 4,
        "phenotype": 5,
        "alias": 6,
        "HPO": 7,
        "tags": 8,
    }

    def __init__(
        self, attributes: list[str] | dict[str, str | int | list[str] | None]
    ) -> None:
        if len(attributes) not in (6, 9):
            # 6 or 9 to allow for ped or ped9 files
            raise ValueError(
                f"Unexpected length of sample attributes list: {attributes}"
            )
        if isinstance(attributes, list):
            self.load_from_list(attributes)
        elif isinstance(attributes, dict):
            self.load_from_dict(attributes)
        else:
            raise ValueError(
                f"Sample attributes should be a list or a dict, got instead:{attributes}"
            )

    def __str__(self):
        return ",".join([str(getattr(self, k)) for k in self.ATTRIBUTES.keys()])

    def load_from_list(self, attributes: list[str]) -> None:
        for k, v in self.ATTRIBUTES.items():
            if v >= 6 and len(attributes) == 6:
                setattr(self, k, [])
            else:
                setattr(self, k, attributes[v])

    def load_from_dict(
        self, attributes: dict[str, str | int | list[str] | None]
    ) -> None:
        """
        Translate from Ped_raw JSON to standard column names
        """
        for k in self.ATTRIBUTES.keys():
            if k == "individual_id":
                setattr(self, k, attributes["id"])
            elif k == "family_id":
                setattr(self, k, attributes["famID"])
            elif k == "paternal_id":
                setattr(self, k, attributes["paternalID"])
            elif k == "maternal_id":
                setattr(self, k, attributes["maternalID"])
            elif k == "HPO":
                setattr(self, k, attributes["HPOList"])
            elif k == "tags":
                setattr(self, k, attributes["starkTags"])
            else:
                setattr(self, k, attributes[k])


class Ped(UserDict):
    """
    UserDict simulates a dictionary. The instance's contents are kept in a regular dictionary, which is accessible via the data attribute.
    See https://docs.python.org/3/library/collections.html#collections.UserDict
    This allows easier use of the Ped object.
        e.g you can directly call
            ped_object["sample_name"]
        instead of
            ped_object.samples["sample_name"] # where samples would have been a dict attribute containing samples
    """

    def __init__(self, ped_file: str | None = None) -> None:
        if not ped_file:
            self.data = {}
            return
        if ped_file.endswith(".json"):
            self.data = self.load_from_json(ped_file)
        elif ped_file.endswith(".ped") or ped_file.endswith(".ped9"):
            self.data = self.load_from_ped9(ped_file)
        else:
            raise ValueError(
                f"Expected file extension: .json, .ped, .ped9 ; got instead: {ped_file}"
            )

    def __str__(self):
        return "\n".join([str(s) for s in self.data.values()])

    def __iter__(self):
        for sample in self.data.values():
            yield sample

    def __contains__(self, value):
        if value in self.data:
            return True
        return False

    def write(self):
        raise NotImplementedError("Only reading is implemented for now")

    def _get_sample(self, sample_name) -> Sample | None:
        try:
            sample = self.data[sample_name]
            return sample
        except KeyError:
            log.warning(f"Sample {sample_name} not found in Ped")
            return None

    def get_maternal_aliases(self, sample_name: str) -> list[str]:
        sample = self._get_sample(sample_name)
        if not sample:
            return []
        parental_aliases = []
        if sample.maternal_id not in ("", None):
            parental_aliases.append(sample.maternal_id)
            if self.data.get(sample.maternal_id, "") not in ("", None):
                parental_aliases += self.data[sample.maternal_id].alias
        return parental_aliases

    def get_paternal_aliases(self, sample_name: str) -> list[str]:
        sample = self._get_sample(sample_name)
        if not sample:
            return []
        parental_aliases = []
        if sample.paternal_id not in ("", None):
            parental_aliases.append(sample.paternal_id)
            if self.data.get(sample.paternal_id, "") not in ("", None):
                parental_aliases += self.data[sample.paternal_id].alias
        return parental_aliases

    def get_samples_from_family(self, family: str) -> list[Sample]:
        return [s for s in self.data.values() if s.family == family]

    def get_non_parental_samples_from_family(self, sample_name: str) -> list[Sample]:
        sample = self._get_sample(sample_name)
        if not sample:
            return []
        family_samples = self.get_samples_from_family(sample.individual_id)
        parental_aliases = self.get_maternal_aliases(
            sample.individual_id
        ) + self.get_paternal_aliases(sample.individual_id)
        return [s for s in family_samples if s.individual_id not in parental_aliases]

    @staticmethod
    def load_from_json(ped_file: str) -> dict[str, Sample]:
        samples = {}
        with open(ped_file, "r") as f:
            ped = json.load(f)
        for sample in ped:
            samples[sample["id"]] = Sample({k: v for k, v in sample.items()})
        return samples

    @staticmethod
    def load_from_ped9(ped_file: str) -> dict[str, Sample]:
        samples = {}
        with open(ped_file, "r") as f:
            for l in f:
                if l.startswith("#"):
                    continue
                if l == "":
                    continue
                l = l.strip("\r\n").split("\t")
                sample = Sample(l)
                samples[sample.individual_id] = sample
        return samples
