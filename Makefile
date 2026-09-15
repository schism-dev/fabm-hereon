# SPDX-FileCopyrightText: 2022-2026 Helmholtz-Zentrum hereon GmbH
# SPDX-License-Identifier: CC0-1.0
# SPDX-FileContributor Carsten Lemmen <carsten.lemmen@hereon.de>
#
# Zenodo metadata workflow for this repository.
#
# This updates the metadata of the EXISTING published record in place
# (same DOI, same version) - it does not create a new version/release.
# See https://developers.zenodo.org/ for background on the two Zenodo
# APIs involved: the classic Deposit API (used here, via zenodraft) and
# the newer native Records API (a different, incompatible schema - not
# handled by this Makefile).
#
# Requires:
#   - zenodraft (https://github.com/zenodraft/zenodraft)
#   - a personal access token from
#       https://zenodo.org/account/settings/applications         (production)
#       https://sandbox.zenodo.org/account/settings/applications (sandbox)
#     exported as ZENODO_ACCESS_TOKEN (or ZENODO_SANDBOX_ACCESS_TOKEN
#     when SANDBOX=1)
#
# Usage:
#   make zenodo-validate                 # offline schema check, no token needed
#   make zenodo-update                   # full flow: unlock, push, publish, verify
#   make zenodo-update SANDBOX=1         # same, against sandbox.zenodo.org
#   make zenodo-verify                   # just print the current live metadata
#
# To publish a NEW version (e.g. backfilling an older git tag as its own
# Zenodo version, with its own source archive):
#   make zenodo-new-version                                  # note the printed version_id
#   make zenodo-archive TAG=v0.2.0                            # builds a zip from that git tag
#   make zenodo-add-file ZENODO_RECORD_ID=<id> FILE=<zip>
#   make zenodo-push ZENODO_RECORD_ID=<id>                    # push .zenodo.json (edit version/
#                                                              # publication_date in it first!)
#   make zenodo-publish ZENODO_RECORD_ID=<id>
#   make zenodo-verify ZENODO_RECORD_ID=<id>
#
# The record/concept ids below are specific to this repository's Zenodo
# deposit (DOI 10.5281/zenodo.21982233, concept DOI 10.5281/zenodo.21982232)
# and only need to change if the record is ever recreated from scratch.

ZENODO_METADATA    ?= .zenodo.json
ZENODO_RECORD_ID   ?= 21982233
ZENODO_CONCEPT_ID  ?= 21982232
ZENODO_ARCHIVE_DIR ?= /tmp/fabm-hereon-archive

ifdef SANDBOX
ZENODO_ZDFLAG := -s
ZENODO_API    := https://sandbox.zenodo.org/api
ZENODO_TOKEN  := $(ZENODO_SANDBOX_ACCESS_TOKEN)
ZENODO_TOKEN_VAR := ZENODO_SANDBOX_ACCESS_TOKEN
else
ZENODO_ZDFLAG :=
ZENODO_API    := https://zenodo.org/api
ZENODO_TOKEN  := $(ZENODO_ACCESS_TOKEN)
ZENODO_TOKEN_VAR := ZENODO_ACCESS_TOKEN
endif

.PHONY: help
help:
	@echo "Targets:"
	@echo "  zenodo-validate      Validate $(ZENODO_METADATA) against the Zenodo schema (no token needed)"
	@echo "  zenodo-unlock        Unlock record $(ZENODO_RECORD_ID) for editing"
	@echo "  zenodo-push          Push $(ZENODO_METADATA) onto the unlocked record"
	@echo "  zenodo-publish       Re-publish the record (finalizes the metadata change)"
	@echo "  zenodo-verify        Print the current live metadata of the record"
	@echo "  zenodo-update        Run unlock, push, publish, verify in sequence"
	@echo "  zenodo-new-version   Create a new draft version under the concept (separate DOI/version)"
	@echo "  zenodo-archive       Build a source zip from a git tag (pass TAG=vX.Y.Z)"
	@echo "  zenodo-add-file      Attach a file to a draft version (pass ZENODO_RECORD_ID= and FILE=)"
	@echo ""
	@echo "Add SANDBOX=1 to any target to run against sandbox.zenodo.org instead."

.PHONY: zenodo-check-token
zenodo-check-token:
	@if [ -z "$(ZENODO_TOKEN)" ]; then \
		echo "error: $(ZENODO_TOKEN_VAR) is not set."; \
		echo "Get a token from $(ZENODO_API:/api=)/account/settings/applications and export it."; \
		exit 1; \
	fi

.PHONY: zenodo-validate
zenodo-validate:
	zenodraft metadata validate $(ZENODO_METADATA)

.PHONY: zenodo-unlock
zenodo-unlock: zenodo-check-token
	@echo "Unlocking record $(ZENODO_RECORD_ID) for editing..."
	@curl -s -X POST -H "Authorization: Bearer $(ZENODO_TOKEN)" \
		"$(ZENODO_API)/deposit/depositions/$(ZENODO_RECORD_ID)/actions/edit" \
		| python3 -m json.tool

.PHONY: zenodo-push
zenodo-push: zenodo-check-token zenodo-validate
	zenodraft metadata update $(ZENODO_ZDFLAG) $(ZENODO_RECORD_ID) $(ZENODO_METADATA)

.PHONY: zenodo-publish
zenodo-publish: zenodo-check-token
	zenodraft deposition publish $(ZENODO_ZDFLAG) $(ZENODO_RECORD_ID)

.PHONY: zenodo-verify
zenodo-verify:
	@echo "Live metadata for record $(ZENODO_RECORD_ID):"
	@curl -s "$(ZENODO_API)/records/$(ZENODO_RECORD_ID)" | python3 -m json.tool

.PHONY: zenodo-update
zenodo-update: zenodo-unlock zenodo-push zenodo-publish zenodo-verify
	@echo "Done. Review the metadata above; the Zenodo web UI may lag a moment behind the API."

.PHONY: zenodo-new-version
zenodo-new-version: zenodo-check-token
	@echo "Creating a new draft version under concept $(ZENODO_CONCEPT_ID)..."
	zenodraft deposition create version $(ZENODO_ZDFLAG) $(ZENODO_CONCEPT_ID)
	@echo "Note the returned version_id above, then typically:"
	@echo "  make zenodo-archive TAG=vX.Y.Z"
	@echo "  make zenodo-add-file ZENODO_RECORD_ID=<version_id> FILE=$(ZENODO_ARCHIVE_DIR)/fabm-hereon-vX.Y.Z.zip"
	@echo "  make zenodo-push ZENODO_RECORD_ID=<version_id>"
	@echo "  make zenodo-publish ZENODO_RECORD_ID=<version_id>"

.PHONY: zenodo-archive
zenodo-archive:
	@if [ -z "$(TAG)" ]; then \
		echo "error: pass TAG=vX.Y.Z (an existing git tag)"; \
		exit 1; \
	fi
	mkdir -p $(ZENODO_ARCHIVE_DIR)
	git archive --format=zip --prefix=fabm-hereon-$(patsubst v%,%,$(TAG))/ $(TAG) \
		-o $(ZENODO_ARCHIVE_DIR)/fabm-hereon-$(TAG).zip
	@echo "Built $(ZENODO_ARCHIVE_DIR)/fabm-hereon-$(TAG).zip from tag $(TAG)"

.PHONY: zenodo-add-file
zenodo-add-file: zenodo-check-token
	@if [ -z "$(FILE)" ]; then \
		echo "error: pass FILE=<path to the archive to attach>"; \
		exit 1; \
	fi
	zenodraft file add $(ZENODO_ZDFLAG) $(ZENODO_RECORD_ID) $(FILE)
