"""
probably a keydb database
key: chr-pos-ref-alt-[all sv fields]
value: json with all annotations
a simple key-value store, no need for complex queries

the "hard" part is to create the initial annotations from howard's databases
ironically, an easy way would be to isolate unknown variants and run them through howard
so this script would more be like a cache for howard

an alternative is to create a dataclasse for each annotation source
    then parse the file into it
    then load it into keydb
there could be a migration script each time the annotation source changes
    that compare the old and new dataclasses and remove/update fields accordingly
maybe those dataclasses could be created from howard's configs?

there could be a master file listing each data source with a hash
    if a data source changes, the annotations are reloaded

if there is a calculation defined in howard, run it after loading the data sources then save the result in the db
also create a dataclass defining the columns created by the calculation

##database structure
keys:
    variant_hash:key_generation [dict of fields:values used to generate the key]
    variant_hash:data_source_1
    variant_hash:data_source_2
    etc...
values: typed key/value pairs of annotations (data will be typed: str, int, float, list, dict)
and also:
    key: data_source_1:version ; values: hash_of_file, last_updated_timestamp, location, scheme (dict based on howard's header)
    key: data_source_2:version
    etc...
    this would store the scheme changes over time which is required for migrations
"""

import json
import redis

class KeyDB:
    """Singleton class to manage KeyDB connection."""
    _instance = None

    def __new__(cls, *args, **kwargs):
        if not cls._instance:
            cls._instance = super(KeyDB, cls).__new__(cls)
        return cls._instance
    
    def __init__(
        self, host: str, port: int = 6379, password: str = "", db: int = 0
    ) -> None:
        if password != "":
            self.conn = redis.Redis(host, port, db, password)
        else:
            self.conn = redis.Redis(host, port, db)

    def set(self, key: str, value: str | dict | list) -> None:
        self.conn.set(key, json.dumps(value))

    def get(self, key: str) -> str | dict | list | None:
        res = self.conn.get(key)
        if res:
            return json.loads(res)
        return None

    def invalidate(self, key: str) -> None:
        raise NotImplementedError()

    def invalidate_all(self) -> None:
        raise NotImplementedError()


def init_db():
    pass

if __name__ == "__main__":
    init_db()