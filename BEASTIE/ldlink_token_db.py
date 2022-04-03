from contextlib import contextmanager
import logging
import os
import sqlite3
import sys
import time


def get_tokens_db(db_path):
    logging.debug(f"Using LDlink token db {db_path}")
    db = sqlite3.connect(db_path, timeout=120)
    db.execute(
        """
CREATE TABLE IF NOT EXISTS ldlink_tokens (
    token TEXT, 
    lock_pid INT DEFAULT NULL,
    lock_time DATE DEFAULT NULL,
    CONSTRAINT pk PRIMARY KEY (token)
)
    """
    )
    return db


def insert_token(db, token):
    print(token)
    db.execute("INSERT OR IGNORE INTO ldlink_tokens VALUES (?, NULL, NULL)", ((token,)))


def try_acquire_token(db):
    with db:
        cur = db.cursor()
        cur.execute(
            """
UPDATE ldlink_tokens
SET lock_time = DATE('now'), lock_pid = ?
WHERE token in (select token from ldlink_tokens WHERE lock_time IS NULL LIMIT 1)
RETURNING token
""",
            ((os.getpid(),)),
        )
        ret = cur.fetchone()
        return ret[0] if ret is not None else None


def release_token(db, token):
    with db:
        db.execute(
            """
UPDATE ldlink_tokens 
SET lock_time = NULL, lock_pid = NULL
WHERE token=?
""",
            ((token,)),
        )


@contextmanager
def acquire_ldlink_token(given_token, db_path, timeout=60 * 60 * 2):
    interval = 10

    if not db_path:
        logging.debug(f"Using passed ldlink token {given_token}")
        try:
            yield given_token
        finally:
            return

    db = get_tokens_db(db_path)
    if given_token:
        logging.debug(f"Inserted ldlink token {given_token} into DB")
        insert_token(db, given_token)

    logging.debug("Begin acquire ldlink_token")
    start = time.time()
    attempts = 0
    while time.time() < start + timeout:
        if attempts % 4 == 0 and attempts != 0:
            logging.debug(
                f"Attempted to acquire token {attempts} times {time.time() - start} seconds"
            )
        attempts += 1
        resource = try_acquire_token(db)
        if resource:
            break
        time.sleep(interval)

    if resource is None:
        logging.error(f"UNABLE TO GET LDLINK TOKEN after {timeout} seconds")
        sys.exit(1)

    try:
        logging.debug(f"Successfully acquired ldlink_token {resource}")
        yield resource
    finally:
        release_token(db, resource)
        logging.debug(f"Released ldlink_token {resource}")
        db.close()
