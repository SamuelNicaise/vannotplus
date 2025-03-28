from collections import UserDict
import json


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

    def is_affected(self) -> bool:
        """
        Check if the sample is affected
        From ped9 specification:
        Affected status should be coded as follows:
        -9 missing
        0 missing
        1 unaffected
        2 affected
        If any value outside of -9,0,1,2 is detected, then the samples are assumed to have phenotype values, interpreted as string phenotype values.
        """
        if self.phenotype not in (None, -9, 0, 1):
            return True
        return False

    def get_parents(self) -> list[str]:
        """
        Returns mother then father ID if both are available

        Returns an empty list if no parents are available
        """
        parents = [self.maternal_id, self.paternal_id]
        parents = [p for p in parents if p not in (None, "")]
        return parents


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

    def get_family_from_sample(self, sample: str) -> str | None:
        if sample not in self.data:
            return None
        return self.data[sample].family_id

    def get_samples_from_family(self, family: str) -> list[Sample]:
        return [s for s in self.data.values() if s.family_id == family]

    def get_children_from_sample(self, sample: str) -> list[Sample]:
        if sample not in self.data:
            raise KeyError(f"Sample {sample} not found in the pedigree.")
        children = []
        for s in self.data.values():
            if s.paternal_id == sample or s.maternal_id == sample:
                children.append(s)
        return children
