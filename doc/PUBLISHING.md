# Publishing the HTML documentation

The existing `index.html` is the home page. GitHub Actions publishes the
contents of `doc/` without a static-site generator. Editor backup files and
macOS desktop metadata are excluded from the deployed artifact.

## First deployment

1. Commit and push the documentation changes and `.github/workflows/pages.yml`
   to `master`.
2. In the GitHub repository, open **Settings > Pages**. Under **Build and
   deployment**, set **Source** to **GitHub Actions**.
3. Open **Actions > Publish documentation > Run workflow**, select `master`,
   and run it. Later pushes changing `doc/` or the workflow deploy automatically.

The site URL is https://juanjosepalacios.github.io/ANT.Gaussian/ .
The ANT.Gaussian section is `ant_g03.html` and its manual is `manual_g03.html`.
The workflow's deployment summary reports the live URL when deployment succeeds.

## Local preview

From the repository root run `python3 -m http.server 8000 --directory doc`,
then open http://localhost:8000/ . Keep asset and navigation URLs relative
so they also work under the GitHub Pages repository path.

GitHub setup reference:
https://docs.github.com/en/pages/getting-started-with-github-pages/using-custom-workflows-with-github-pages
