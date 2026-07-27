OFFICIAL_ENV := $(shell sed -n '1p' envs/official-name.txt)
ENV ?= $(shell ./run_in_env.sh python scripts/envctl.py name 2>/dev/null | tail -n 1)

.PHONY: env-bootstrap env-render env-lock env-pip-lock env-create env-install env-check env-test env-snapshot

env-bootstrap:
	./run_in_env.sh --create-env --env $(OFFICIAL_ENV) --file envs/locks/hic-basic-linux-64.lock
	./run_in_env.sh --env $(OFFICIAL_ENV) --network internet python scripts/envctl.py install-editable
	./run_in_env.sh --env $(OFFICIAL_ENV) python scripts/envctl.py verify

env-render:
	./run_in_env.sh python scripts/envctl.py render

env-lock: env-render
	./run_in_env.sh --network internet python scripts/envctl.py lock

env-pip-lock: env-render
	./run_in_env.sh --network internet python scripts/envctl.py lock-pip

env-create:
	./run_in_env.sh python scripts/envctl.py create --env $(ENV)

env-install:
	./run_in_env.sh --env $(ENV) --network internet python scripts/envctl.py install-editable

env-check:
	./run_in_env.sh --env $(ENV) python scripts/envctl.py verify

env-test:
	./run_in_env.sh --env $(ENV) pytest -q tests/test_envctl.py tests/test_api_cli.py tests/test_state_track.py tests/test_interpolate.py

env-snapshot:
	./run_in_env.sh --env $(ENV) python scripts/envctl.py snapshot --output .env-backups/$(ENV).json
