# Docker runtime

This `Dockerfile` demonstrates how to create a simple sawfish runtime docker image.

Because sawfish is a single static binary with no dependencies, it is not distributed as a docker image and
it is not recommended to run it through docker for typical workflows. This Docker template is intended for
debug and test scenarios only.
