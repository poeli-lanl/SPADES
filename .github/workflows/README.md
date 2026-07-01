# GitHub Actions Workflows

This directory contains CI/CD workflows for the SPADES project.

## Workflows

### 1. `test.yml` - Run Tests
**Triggers:** Push to main/master/develop, Pull Requests, Manual dispatch

**Purpose:** Runs the test suite to validate functionality.

**Steps:**
- Sets up Conda environment from `environment.yml`
- Runs `post_install.sh` to configure the environment
- Executes `scripts/check_install.sh` to verify installation
- Runs test suite if available
- Uploads test results as artifacts

### 2. `docker-test.yml` - Docker Build Test
**Triggers:** Push/PR to main/master/develop (when Dockerfile or related files change), Manual dispatch

**Purpose:** Tests Docker image building without publishing.

**Steps:**
- Builds Docker image using Buildx
- Runs `check_install.sh` inside the container
- Reports image size

### 3. `docker-build.yml` - Build and Push Docker Image
**Triggers:** New release published, Manual dispatch

**Purpose:** Automatically builds and publishes Docker images to Docker Hub when a new release is created.

**Steps:**
- Builds multi-architecture Docker image (AMD64)
- Tags with version from release tag and `latest`
- Pushes to Docker Hub (`poeli/spades-g2`)
- Tests the built image
- Updates Docker Hub repository description

## Setup Instructions

### Required GitHub Secrets

To enable Docker Hub publishing, you need to configure the following secrets in your GitHub repository:

1. Go to your repository on GitHub
2. Navigate to **Settings** → **Secrets and variables** → **Actions**
3. Click **New repository secret** and add:

#### `DOCKERHUB_USERNAME`
- Your Docker Hub username (e.g., `poeli`)

#### `DOCKERHUB_TOKEN`
- A Docker Hub access token (recommended) or password
- To create an access token:
  1. Log in to [Docker Hub](https://hub.docker.com)
  2. Go to **Account Settings** → **Security** → **Access Tokens**
  3. Click **New Access Token**
  4. Give it a description (e.g., "GitHub Actions SPADES")
  5. Set permissions to **Read & Write**
  6. Copy the token and add it to GitHub secrets

### Creating a Release

To trigger the Docker build workflow:

1. **Via GitHub UI:**
   - Go to your repository's **Releases** page
   - Click **Draft a new release**
   - Create a new tag (e.g., `v1.2.1` or `1.2.1`)
   - Fill in release title and description
   - Click **Publish release**

2. **Via Git Command Line:**
   ```bash
   git tag -a v1.2.1 -m "Release version 1.2.1"
   git push origin v1.2.1
   ```
   Then create the release on GitHub using the pushed tag.

3. **Manual Trigger:**
   - Go to **Actions** → **Build and Push Docker Image**
   - Click **Run workflow**
   - Enter a custom tag name
   - Click **Run workflow**

### Docker Image Tags

When a release is published:
- `poeli/spades-g2:X.Y.Z` (version from the release tag)
- `poeli/spades-g2:latest` (always points to the latest release)

### Testing Locally

Before pushing, you can test the workflows locally using [act](https://github.com/nektos/act):

```bash
# Install act
brew install act  # macOS
# or
curl https://raw.githubusercontent.com/nektos/act/master/install.sh | sudo bash  # Linux

# Test the test workflow
act push -W .github/workflows/test.yml

# Test Docker build (without pushing)
act push -W .github/workflows/docker-test.yml
```

## Workflow Status

You can view workflow runs and their status at:
`https://github.com/[your-username]/[your-repo]/actions`

## Troubleshooting

### Docker Build Fails
- Check Dockerfile syntax
- Verify `environment.yml` is valid
- Ensure all required files are present
- Check workflow logs for specific errors

### Docker Push Fails
- Verify Docker Hub credentials are correct
- Check that `DOCKERHUB_TOKEN` has read/write permissions
- Ensure the Docker Hub repository exists

### Tests Fail
- Check test database availability
- Verify Conda environment setup
- Review test logs in the Actions tab
- Download test artifacts for debugging

## Contributing

When modifying workflows:
1. Test changes in a feature branch
2. Use `workflow_dispatch` for manual testing
3. Review logs carefully
4. Update this README if adding new workflows
