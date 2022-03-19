import sqlite3
from time import sleep


def try_acquire_token(db):
    with db:
        cur = db.cursor()
        cur.execute(
            """
UPDATE ldlink_tokens
SET lock_time = DATE('now')
WHERE lock_time IS NULL
RETURNING token
LIMIT 1
"""
        )
        ret = cur.fetchone()
        return ret[0] if ret is not None else None


def release_token(db, token):
    with db:
        db.execute(
            """
UPDATE ldlink_tokens 
SET lock_time = NULL 
WHERE token=?
""",
            token,
        )


db = sqlite3.connect("test_ldlink_tokens.db")
db.execute(
    """
CREATE TABLE IF NOT EXISTS ldlink_tokens (
    token TEXT, 
    lock_time DATE DEFAULT NULL,
    CONSTRAINT pk PRIMARY KEY (token)
)
"""
)

for i in range(5):
    token = try_acquire_token(db)
    if token is not None:
        print("got token", token)
        sleep(3)
        release_token(db, token)
        break
    else:
        print("unable to get lock, waiting to try again")
        sleep(1)

db.close()
