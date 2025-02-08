# CI Steps
#
# Considerations
# - Easy to debug: show the command being run
# - Leverage CI features: Only run individual steps so we can use features like reporting elapsed time per step

ARGS?=--workspace
TOOLCHAIN_TARGET ?=
ifneq (${TOOLCHAIN_TARGET},)
  ARGS+=--target ${TOOLCHAIN_TARGET}
endif

STABLE?=1.84

_FEATURES = minimal default debug release
_FEATURES_minimal = --no-default-features --features "std"
_FEATURES_default =
_FEATURES_debug = ${_FEATURES_default} --features debug
_FEATURES_release = ${_FEATURES_default} --release

check-%:
	cargo check ${_FEATURES_${@:check-%=%}} --all-targets ${ARGS}

build-%:
	cargo test ${_FEATURES_${@:build-%=%}} --all-targets --no-run ${ARGS}

test-%:
	cargo test ${_FEATURES_${@:test-%=%}} ${ARGS}

clippy-%:
	cargo clippy ${_FEATURES_${@:clippy-%=%}} ${ARGS} --all-targets -- -D warnings -A deprecated

doc:
	cargo doc --workspace --all-features --no-deps --document-private-items