BUILD    := $(shell git describe --tags --always --dirty)
DIR      := $(shell pwd)

REGISTRY := seglh
APP      := alu_analysis
IMG      := $(REGISTRY)/$(APP)
IMG_VERSIONED := $(IMG):$(BUILD)
IMG_LATEST    := $(IMG):latest
TAR      := $(DIR)/$(REGISTRY)-$(APP)-$(BUILD).tar.gz

.PHONY: all version build push cleanbuild

all: build

version:
	@echo "Build version: $(BUILD)"

build: version
	docker build -f Dockerfile \
		-t $(IMG_VERSIONED) .
	docker tag $(IMG_VERSIONED) $(IMG_LATEST)
	docker save $(IMG_VERSIONED) | gzip > $(TAR)
	@echo "Saved $(TAR)"

push:
	docker push $(IMG_VERSIONED)
	docker push $(IMG_LATEST)

cleanbuild: version
	docker build -f Dockerfile \
		--no-cache \
		-t $(IMG_VERSIONED) .
	docker tag $(IMG_VERSIONED) $(IMG_LATEST)
	docker save $(IMG_VERSIONED) | gzip > $(TAR)
	@echo "Saved $(TAR)"