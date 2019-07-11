pipeline {
  agent any
  stages {
    stage('Docker build image') {
      steps {
        sh '''docker build -t mdstudio/smartcyp .'''
      }
    }
  }
}