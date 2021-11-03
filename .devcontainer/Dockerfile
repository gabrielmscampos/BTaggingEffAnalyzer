ARG VARIANT="3.10-bullseye"
FROM mcr.microsoft.com/vscode/devcontainers/python:0-${VARIANT}

# Update pip
RUN pip install --no-cache-dir -U pip

# Install poetry
ARG POETRY_VERSION="none"
USER vscode
RUN curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/install-poetry.py | python - --version=${POETRY_VERSION}
ENV PATH "/home/vscode/.poetry/bin:$PATH"
RUN /home/vscode/.local/bin/poetry config virtualenvs.create false \
    && sudo chown -R vscode /usr/local/src

# Return to root user after installing Poetry
USER root

# [Choice] Node.js version: none, lts/*, 16, 14, 12, 10
ARG NODE_VERSION="none"
RUN if [ "${NODE_VERSION}" != "none" ]; then echo "Installing Node.js..." \
    && vscode -c "umask 0002 && . /usr/local/share/nvm/nvm.sh && nvm install ${NODE_VERSION} 2>&1"; fi

# [Choice] Install docker 
ARG DOCKER="no"
RUN if [ "${DOCKER}" != "no" ]; then echo "Installing Docker engine..." \
    && curl -fsSL https://get.docker.com -o get-docker.sh \
    && apt-get install -y wget \
    && sh get-docker.sh \
    && usermod -aG docker root \
    && apt-get autoremove; fi

# [Choice] Install docker compose
ARG DOCKER_COMPOSE="no"
RUN if [ "${DOCKER_COMPOSE}" != "no" ]; then echo "Installing Docker Compose..." \
    && LATEST_COMPOSE_VERSION=$(curl -sSL "https://api.github.com/repos/docker/compose/releases/latest" | grep -o -P '(?<="tag_name": ").+(?=")') \
    && curl -sSL "https://github.com/docker/compose/releases/download/${LATEST_COMPOSE_VERSION}/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose \
    && chmod +x /usr/local/bin/docker-compose; fi

# [Option] Install aws cli
ARG AWS_CLI="no"
RUN if [ "${AWS_CLI}" != "no" ]; then echo "Installing aws cli..." \
    && curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" \
    && unzip awscliv2.zip \
    && sudo ./aws/install \
    && rm -rf ./awscliv2.zip ./aws; fi